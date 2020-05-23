//
// Seidr - Create and operate on gene crowd networks
// Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
//
// This file is part of Seidr.
//
// Seidr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Seidr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
//

#include <BSlogger.hpp>
#include <Serialize.h>
#include <common.h>
#include <fs.h>
#include <parallel_control.h>
#include <threshold.h>

#undef DEBUG

#if defined(SEIDR_PSTL)
#include <pstl/algorithm>
#include <pstl/execution>
#include <tbb/task_scheduler_init.h>
#else
#include <algorithm>
#endif
#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <networkit/global/ClusteringCoefficient.h>
#include <networkit/graph/Graph.h>
#include <string>
#include <vector>

namespace po = boost::program_options;

bool
reverse_sort(const double& a, const double& b)
{
  return a > b;
}

struct thresh_t
{
  std::map<std::string, double> data = {
    { "t", 0 },  { "ne", 0 },  { "ner", 0 }, { "nn", 0 }, { "nnr", 0 },
    { "sf", 0 }, { "sfr", 0 }, { "cc", 0 },  { "ccr", 0 }
  };
};

bool
ne_sort(const thresh_t& a, const thresh_t& b)
{
  return boost::numeric_cast<uint32_t>(a.data.at("ne")) >
         boost::numeric_cast<uint32_t>(b.data.at("ne"));
}
bool
nn_sort(const thresh_t& a, const thresh_t& b)
{
  return boost::numeric_cast<uint32_t>(a.data.at("nn")) >
         boost::numeric_cast<uint32_t>(b.data.at("nn"));
}
bool
sf_sort(const thresh_t& a, const thresh_t& b)
{
  return a.data.at("sf") > b.data.at("sf");
}
bool
cc_sort(const thresh_t& a, const thresh_t& b)
{
  // NOLINTNEXTLINE: "Magic number" required for math
  return std::abs(0.5 - a.data.at("cc")) < std::abs(0.5 - b.data.at("cc"));
}
bool
t_sort(const thresh_t& a, const thresh_t& b)
{
  return a.data.at("t") > b.data.at("t");
}

void
get_ranks(std::vector<thresh_t>& v, const std::string& s)
{
  std::string tr;
  tr = s + "r";

  if (s == "ne") {
    SORTWCOMP(v.begin(), v.end(), ne_sort);
  } else if (s == "nn") {
    SORTWCOMP(v.begin(), v.end(), nn_sort);
  } else if (s == "sf") {
    SORTWCOMP(v.begin(), v.end(), sf_sort);
  } else if (s == "cc") {
    SORTWCOMP(v.begin(), v.end(), cc_sort);
  }

  double rank = 0;
  v[0].data[tr] = rank;

  for (uint64_t cnt = 1; cnt < v.size(); cnt++) {
    if (almost_equal(v[cnt].data[s], v[cnt - 1].data[s])) {
      v[cnt].data[tr] = rank;
    } else {
      rank++;
      v[cnt].data[tr] = rank;
    }
  }
}

thresh_t
check_one(std::vector<MiniEdge>& ev,
          double thresh,
          uint64_t& stop,
          bool rank,
          std::ostream& out)
{
  logger log(std::cerr, "threshold");

  thresh_t ret;

  log(LOG_DEBUG) << "Processing threshold " << thresh << '\n';
  if (rank) {
    while (stop < ev.size() && ev[stop].s < thresh) {
      stop++;
    }
  } else {
    while (stop < ev.size() && ev[stop].s > thresh) {
      stop++;
    }
  }
  stop = stop > 0 ? stop - 1 : 0;

  std::map<uint64_t, uint64_t> id_remap;

  log(LOG_DEBUG) << "Remapping IDs\n";
  for (uint64_t i = 0; i < stop; i++) {
    id_remap[ev[i].i] = 0;
    id_remap[ev[i].j] = 0;
  }

  uint64_t nid = 0;
  for (auto& x : id_remap) {
    x.second = nid++;
  }

  log(LOG_DEBUG) << "Creating networkit graph\n";

  NetworKit::Graph g(id_remap.size(), false, false);
  for (uint64_t i = 0; i < stop; i++) {
    g.addEdge(id_remap[ev[i].i], id_remap[ev[i].j]);
  }

  auto ne = g.numberOfEdges();
  auto nv = g.numberOfNodes();

  log(LOG_DEBUG) << "Calculating stats\n";

  std::map<int64_t, uint64_t> deg_map;
  for (uint64_t i = 0; i < nv; i++) {
    double x = g.degree(i);
    auto it = deg_map.find(x);
    if (it == deg_map.end()) {
      deg_map[x] = 0;
    }
    deg_map[x]++;
  }

  std::vector<double> x;
  std::vector<double> y;
  double dnv = nv;
  for (auto& it : deg_map) {
    double d = log10(boost::numeric_cast<double>(it.first + 1));
    double pd = log10(boost::numeric_cast<double>(it.second + 1) /
                      boost::numeric_cast<double>(dnv + 1));
    x.push_back(d);
    y.push_back(pd);
  }

  arma::mat c = arma::cor(arma::conv_to<arma::vec>::from(x),
                          arma::conv_to<arma::vec>::from(y));
  double rsqu = arma::as_scalar(arma::pow(c, 2));

  double trans = NetworKit::ClusteringCoefficient::exactGlobal(g);

  out << thresh << '\t' << nv << '\t' << ne << '\t' << rsqu << '\t' << trans
      << '\n';

  ret.data["t"] = thresh;
  ret.data["ne"] = ne;
  ret.data["nn"] = nv;
  ret.data["sf"] = rsqu;
  ret.data["cc"] = trans;
  return ret;
}

std::vector<double>
make_steps(double min, double max, uint32_t nsteps, bool rank)
{
  std::vector<double> ret;
  double x = nsteps;
  double alpha = (max - min) / x;
  double beta = min;
  for (uint32_t i = 0; i < nsteps; i++) {
    ret.push_back(beta);
    beta += alpha;
  }
  ret.push_back(max);
  if (rank) {
    SORTWCOMP(ret.begin(), ret.end(), reverse_sort);
  }
  return ret;
}

int
threshold(const std::vector<std::string>& args)
{
  logger log(std::cerr, "threshold");

  try {

    seidr_threshold_param_t param;

    po::options_description umbrella(
      "Pick hard network threshold according to topology");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("-"),
      "Output file name ['-' for stdout]");

    po::options_description ompopt("OpenMP Options");
    ompopt.add_options()(
      "threads,O",
      po::value<int>(&param.nthreads)->default_value(GET_MAX_PSTL_THREADS()),
      "Number of OpenMP threads for parallel sorting");

    po::options_description topt("Threshold Options");
    topt.add_options()("min,m",
                       po::value<double>(&param.min)->default_value(0.0, "0"),
                       "Lowest threshold value to check")(
      "max,M",
      po::value<double>(&param.max)->default_value(1.0, "0"),
      "Highest threshold value to check")(
      "index,i",
      po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
      "Score column to use as edge weights")(
      "nsteps,n",
      po::value<uint32_t>(&param.nsteps)
        ->default_value(SEIDR_THRESHOLD_DEF_NSTEP),
      "Number of breaks to create for testing");

    po::options_description fopt("Formatting Options");
    fopt.add_options()("precision,p",
                       po::value<uint16_t>(&param.precision)
                         ->default_value(SEIDR_THRESHOLD_DEF_PREC),
                       "Number of decimal points to print");

    po::options_description req("Required [can be positional]");
    req.add_options()("in-file",
                      po::value<std::string>(&param.infile)->required(),
                      "Input SeidrFile");

    umbrella.add(req).add(topt).add(fopt).add(ompopt).add(opt);

    po::positional_options_description p;
    p.add("in-file", 1);

    po::variables_map vm;
    po::store(
      po::command_line_parser(args).options(umbrella).positional(p).run(), vm);

    if (vm.count("help") != 0) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);

    set_pstl_threads(param.nthreads);
    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);

    if (param.outfile != "-") {
      param.outfile = to_absolute(param.outfile);
      if (!param.force) {
        assert_no_overwrite(param.outfile);
      }
      assert_dir_is_writeable(dirname(param.outfile));
    }

    std::shared_ptr<std::ostream> out;
    if (param.outfile == "-") {
      out = std::shared_ptr<std::ostream>(&std::cout, no_delete);
    } else {
      out =
        std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    SeidrFile rf(param.infile.c_str());
    rf.open("r");

    SeidrFileHeader h;
    h.unserialize(rf);

    make_tpos(param.tpos, h);

    log(LOG_INFO) << "Starting analysis\n";
    log(LOG_INFO) << "Creating vector of thresholds\n";
    std::vector<double> steps =
      make_steps(param.min, param.max, param.nsteps, param.trank);

    log(LOG_INFO) << "Reading network with "
                  << (param.trank ? "ranks < " : "scores > ") << steps[0]
                  << '\n';
    std::vector<MiniEdge> edges =
      read_network_minimal(h, rf, steps[0], param.tpos, param.trank);
    log(LOG_INFO) << "Read " << edges.size() << " edges\n";
    log(LOG_INFO) << "Sorting network\n";

    if (param.trank) {
      SORTWCOMP(edges.begin(), edges.end(), me_rank_sort);
    } else {
      SORTWCOMP(edges.begin(), edges.end(), me_score_sort);
    }

    uint64_t stop = 0;

    std::vector<thresh_t> tvec;

    log(LOG_INFO) << "Starting topology assesment\n";
    out->precision(param.precision);
    for (auto it = steps.rbegin(); it != steps.rend(); it++) {
      log(LOG_INFO) << "Current threshold:" << *it << "\n";
      tvec.push_back(check_one(edges, *it, stop, param.trank, *out));
    }

    get_ranks(tvec, "nn");
    get_ranks(tvec, "ne");
    get_ranks(tvec, "sf");
    get_ranks(tvec, "cc");

    double cut =
      tvec[0].data["nnr"] + tvec[0].data["sfr"] + tvec[0].data["ccr"];
    thresh_t final = tvec[0];
    for (auto& t : tvec) {
      if (t.data["nnr"] + t.data["sfr"] + t.data["ccr"] < cut) {
        final = t;
      }
    }

    log(LOG_INFO) << "Suggested threshold: " << final.data["t"] << "("
                  << final.data["nnr"] << ',' << final.data["sfr"] << ','
                  << final.data["ccr"] << ")" << '\n';

    rf.close();

  } catch (const po::error& e) {
    log(LOG_ERR) << "[Argument Error]: " << e.what() << '\n';
    return 1;
  } catch (const std::runtime_error& e) {
    log(LOG_ERR) << "[Runtime Error]: " << e.what() << '\n';
    return 1;
  } catch (const std::exception& e) {
    log(LOG_ERR) << "[Generic Error]: " << e.what() << '\n';
    return 1;
  }

  return 0;
}
