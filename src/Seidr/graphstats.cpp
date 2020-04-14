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
#include <graphstats.h>
#include <stats.h>

#include <algorithm>
#include <armadillo>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace po = boost::program_options;

#undef DEBUG

#include <networkit/base/Algorithm.h>
#include <networkit/components/ConnectedComponents.h>
#include <networkit/distance/APSP.h>
#include <networkit/global/ClusteringCoefficient.h>
#include <networkit/graph/Graph.h>

std::pair<double, double>
SFT(const NetworKit::Graph& g, uint64_t n_bins = 200)
{
  uint64_t nn = g.numberOfNodes();
  double nnd = boost::numeric_cast<double>(nn);

  arma::vec deg_vec(nn, arma::fill::zeros);

  for (uint64_t i = 0; i < nn; i++) {
    deg_vec(i) = g.degree(i);
  }

  double deg_sum = arma::as_scalar(arma::accu(deg_vec));
  double avg_deg = deg_sum / nnd;

  arma::vec centers =
    arma::linspace<arma::vec>(arma::min(deg_vec), arma::max(deg_vec), n_bins);

  arma::vec hist = arma::conv_to<arma::vec>::from(arma::hist(deg_vec, centers));

  centers = arma::log10(centers);
  hist = arma::log10(hist + 1); // Add 1 to avoid log10(0)

  // For the most part R^2 of a linear model is equal to the correlation squared
  double rsq = arma::as_scalar(arma::pow(arma::cor(centers, hist), 2));

  return std::pair<double, double>(avg_deg, rsq);
}

double
avg_w_deg(const NetworKit::Graph& g)
{

  uint64_t nn = g.numberOfNodes();
  double nnd = boost::numeric_cast<double>(nn);

  arma::vec deg_vec(nn, arma::fill::zeros);

  for (uint64_t i = 0; i < nn; i++) {
    deg_vec(i) = g.weightedDegree(i);
  }

  double deg_sum = arma::as_scalar(arma::accu(deg_vec));
  double avg_deg = deg_sum / nnd;

  return avg_deg;
}

std::pair<double, double>
diameter(const std::vector<std::vector<double>>& sps)
{
  double max = 0;
  double np = 0;
  double sum = 0;

  for (const auto& node : sps) {
    for (const auto& pathlength : node) {
      if (!almost_equal(pathlength, std::numeric_limits<double>::max()) &&
          std::isfinite(pathlength)) {
        if (pathlength > max) {
          max = pathlength;
        }
        sum += pathlength;
        np++;
      }
    }
  }

  return std::pair<double, double>(max, sum / np);
}

int
graphstats(int argc, char** argv)
{
  logger log(std::cerr, "graphstats");

  const char* args[argc - 1];
  std::string pr(argv[0]);
  pr += " graphstats";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++)
    args[i - 1] = argv[i];
  argc--;

  seidr_graphstat_param_t param;

  po::options_description umbrella("Calculate graph level network statistics");

  po::options_description opt("Common Options");
  opt.add_options()("force,f",
                    po::bool_switch(&param.force)->default_value(false),
                    "Force overwrite output file if it exists")(
    "help,h", "Show this help message")(
    "outfile,o",
    po::value<std::string>(&param.outfile)->default_value("-"),
    "Output file name ['-' for stdout]");

  po::options_description gsopt("Graphstat Options");
  gsopt.add_options()(
    "index,i",
    po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
    "Index of scores that should be used as weights");

  po::options_description req("Required Options [can be positional]");
  req.add_options()("in-file",
                    po::value<std::string>(&param.infile)->required(),
                    "Input SeidrFile");

  umbrella.add(req).add(gsopt).add(opt);

  po::positional_options_description p;
  p.add("in-file", 1);

  po::variables_map vm;
  po::store(
    po::command_line_parser(argc, args).options(umbrella).positional(p).run(),
    vm);

  if (vm.count("help") || argc == 1) {
    std::cerr << umbrella << '\n';
    return 22;
  }

  try {
    po::notify(vm);
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
  } catch (std::runtime_error& except) {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  std::shared_ptr<std::ostream> out;
  if (param.outfile == "-") {
    out.reset(&std::cout, [](...) {});
  } else {
    out.reset(new std::ofstream(param.outfile.c_str()));
  }

  SeidrFile rf(param.infile.c_str());
  rf.open("r");

  SeidrFileHeader h;
  h.unserialize(rf);

  make_tpos(param.tpos, h);

  log(LOG_INFO) << "Reading graph\n";
  std::vector<SeidrFileEdge> ev = read_network(
    h, rf, -std::numeric_limits<double>::infinity(), param.tpos, false);

  log(LOG_INFO) << "Creating Networkit graph\n";
  NetworKit::Graph g(h.attr.nodes, true, false);
  for (auto& e : ev)
    g.addEdge(e.index.i, e.index.j, e.scores[param.tpos].s);

  (*out) << "Number of Nodes:\t" << g.numberOfNodes() << '\n';
  (*out) << "Number of Edges:\t" << g.numberOfEdges() << '\n';

  log(LOG_INFO) << "Calculating connected components\n";
  NetworKit::ConnectedComponents ccn(g);
  ccn.run();
  uint64_t comp = ccn.numberOfComponents();
  (*out) << "Number of Connected Components:\t" << comp << '\n';

  log(LOG_INFO) << "Calculating clustering coefficient\n";
  NetworKit::ClusteringCoefficient ccoef;
  double gccf = ccoef.exactGlobal(g);
  (*out) << "Global clustering coefficient:\t" << gccf << '\n';

  log(LOG_INFO) << "Calculating scale-freeness\n";
  auto sft = SFT(g);
  (*out) << "Scale free fit:\t" << sft.second << '\n';
  (*out) << "Average degree:\t" << sft.first << '\n';
  (*out) << "Average weighted degree:\t" << avg_w_deg(g) << '\n';

  log(LOG_INFO) << "Calculating network diameter\n";
  NetworKit::APSP apsp(g);
  apsp.run();
  auto all_dist = apsp.getDistances();
  auto dia = diameter(all_dist);
  (*out) << "Network diameter:\t" << dia.first << '\n';
  (*out) << "Average path length:\t" << dia.second << '\n';

  return 0;
}