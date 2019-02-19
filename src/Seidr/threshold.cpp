#include <common.h>
#include <threshold.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>

#undef DEBUG

#include <networkit/graph/Graph.h>
#include <networkit/global/ClusteringCoefficient.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <tclap/CmdLine.h>
#include <cmath>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>
#include <boost/numeric/conversion/cast.hpp>

bool reverse_sort (const double& a, const double& b) {return a > b;}

struct thresh_t {
  std::map<std::string, double> data = {{"t", 0}, {"ne", 0}, {"ner", 0},
    {"nn", 0}, {"nnr", 0}, {"sf", 0},
    {"sfr", 0}, {"cc", 0}, {"ccr", 0}
  };
};

bool ne_sort(const thresh_t& a, const thresh_t& b) {
  return
    boost::numeric_cast<uint32_t>(a.data.at("ne")) >
    boost::numeric_cast<uint32_t>(b.data.at("ne"));
}
bool nn_sort(const thresh_t& a, const thresh_t& b) {
  return
    boost::numeric_cast<uint32_t>(a.data.at("nn")) >
    boost::numeric_cast<uint32_t>(b.data.at("nn"));
}
bool sf_sort(const thresh_t& a, const thresh_t& b) {
  return
    a.data.at("sf") > b.data.at("sf");
}
bool cc_sort(const thresh_t& a, const thresh_t& b) {
  return
    std::abs(0.5 - a.data.at("cc")) <
    std::abs(0.5 - b.data.at("cc"));
}
bool t_sort(const thresh_t& a, const thresh_t& b) {
  return
    a.data.at("t") > b.data.at("t");
}

void get_ranks(std::vector<thresh_t>& v, std::string s)
{
  std::string tr;
  tr = s + "r";

  if (s == "ne")
    std::sort(v.begin(), v.end(), ne_sort);
  else if (s == "nn")
    std::sort(v.begin(), v.end(), nn_sort);
  else if (s == "sf")
    std::sort(v.begin(), v.end(), sf_sort);
  else if (s == "cc")
    std::sort(v.begin(), v.end(), cc_sort);


  double rank = 0;
  v[0].data[tr] = rank;

  for (uint64_t cnt = 1; cnt < v.size(); cnt++)
  {
    if (almost_equal(v[cnt].data[s], v[cnt - 1].data[s]))
    {
      v[cnt].data[tr] = rank;
    }
    else
    {
      rank++;
      v[cnt].data[tr] = rank;
    }
  }

}

thresh_t check_one(std::vector<MiniEdge>& ev, double thresh, uint64_t& stop,
                   bool rank)
{
  logger log(std::cerr, "threshold");

  thresh_t ret;

  log(LOG_DEBUG) << "Processing threshold " << thresh << '\n';
  if (rank)
  {
    while (stop < ev.size() && ev[stop].s < thresh)
    {
      stop++;
    }
  }
  else
  {
    while (stop < ev.size() && ev[stop].s > thresh)
    {
      stop++;
    }
  }
  stop = stop > 0 ? stop - 1 : 0;


  std::map<uint64_t, uint64_t> id_remap;

  log(LOG_DEBUG) << "Remapping IDs\n";
  for (uint64_t i = 0; i < stop; i++)
  {
    id_remap[ev[i].i] = 0;
    id_remap[ev[i].j] = 0;
  }

  uint64_t nid = 0;
  for (auto& x : id_remap)
    x.second = nid++;

  log(LOG_DEBUG) << "Creating networkit graph\n";

  NetworKit::Graph g(id_remap.size(), false, false);
  for (uint64_t i = 0; i < stop; i++)
    g.addEdge(id_remap[ev[i].i], id_remap[ev[i].j]);

  auto ne = g.numberOfEdges();
  auto nv = g.numberOfNodes();

  log(LOG_DEBUG) << "Calculating stats\n";

  std::map<long int, uint64_t> deg_map;
  for (uint64_t i = 0; i < nv; i++)
  {
    double x = g.degree(i);
    auto it = deg_map.find(x);
    if (it == deg_map.end())
    {
      deg_map[x] = 0;
    }
    deg_map[x]++;
  }

  std::vector<double> x, y;
  double dnv = nv;
  for (auto& it : deg_map)
  {
    double d = log10(it.first + 1);
    double pd = log10((it.second + 1) / (dnv + 1));
    x.push_back(d);
    y.push_back(pd);
  }

  double c0 = 0, c1 = 0, cov00 = 0, cov01 = 0, cov11 = 0, sumsq = 0, tss = 0, rsqu = 0;

  gsl_fit_linear (&y[0], 1, &x[0], 1, x.size(),
                  &c0, &c1, &cov00, &cov01, &cov11,
                  &sumsq);
  tss = gsl_stats_tss(&x[0], 1, x.size());
  rsqu = 1 - sumsq / tss;

  NetworKit::ClusteringCoefficient cc;
  double trans = cc.exactGlobal(g);

  std::cout << thresh << '\t' << nv << '\t'
            << ne << '\t' << rsqu << '\t' << trans << '\n';

  ret.data["t"] = thresh;
  ret.data["ne"] = ne;
  ret.data["nn"] = nv;
  ret.data["sf"] = rsqu;
  ret.data["cc"] = trans;
  return ret;
}

std::vector<double> make_steps(double min, double max, uint32_t nsteps,
                               bool rank)
{
  std::vector<double> ret;
  double x = nsteps;
  double alpha = (max - min) / x;
  double beta = min;
  for (uint32_t i  = 0; i < nsteps; i++)
  {
    ret.push_back(beta);
    beta += alpha;
  }
  ret.push_back(max);
  if (rank)
  {
    std::sort(ret.begin(), ret.end(), reverse_sort);
  }
  return ret;
}

int threshold(int argc, char ** argv)
{
  logger log(std::cerr, "threshold");
  double min;
  double max;
  uint32_t nsteps;
  uint32_t tpos;
  std::string infile;
  bool trank;
  uint16_t precision;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " threshold";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Pick network threshold according to topology", ' ',
                       version);

    TCLAP::ValueArg<double>
    arg_min("m", "min", "Lowest edge weight to check", false,
            0, "0");
    cmd.add(arg_min);

    TCLAP::ValueArg<double>
    arg_max("M", "max", "Highest edge weight to check", false,
            1, "1");
    cmd.add(arg_max);

    TCLAP::ValueArg<uint32_t>
    arg_index("i", "index", "Apply threshold on this index", false,
              0, "last column");
    cmd.add(arg_index);

    TCLAP::ValueArg<uint32_t>
    arg_nsteps("n", "nsteps", "Number of steps to test", false,
               100, "100");
    cmd.add(arg_nsteps);

    TCLAP::ValueArg<uint16_t>
    arg_precision("p", "precision", "Number of digits after comma", false,
                  8, "8");
    cmd.add(arg_precision);

    TCLAP::SwitchArg
    arg_r("r", "threshold-rank",
          "Apply threshold value to rank rather than score", cmd, false);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true,
               "", "");

    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    min = arg_min.getValue();
    max = arg_max.getValue();
    nsteps = arg_nsteps.getValue();
    tpos = arg_index.getValue();
    trank = arg_r.getValue();
    precision = arg_precision.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return EINVAL;
  }

  try
  {
    infile = to_absolute(infile);
    if (! file_can_read(infile))
    {
      throw std::runtime_error("Can't read file: " + infile);
    }
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

  SeidrFile rf(infile.c_str());
  rf.open("r");

  SeidrFileHeader h;
  h.unserialize(rf);

  if (tpos == 0)
    tpos = h.attr.nalgs - 1;
  else
    tpos--;

  log(LOG_INFO) << "Starting analysis\n";
  log(LOG_INFO) << "Creating vector of thresholds\n";
  std::vector<double> steps = make_steps(min, max, nsteps, trank);

  log(LOG_INFO) << "Reading network with "
                << (trank ? "ranks < " : "scores > ")
                << steps[0] << '\n';
  std::vector<MiniEdge> edges = read_network_minimal(h, rf, steps[0],
                                tpos, trank);
  log(LOG_INFO) << "Read " << edges.size() << " edges\n";
  log(LOG_INFO) << "Sorting network\n";

  if (trank)
  {
    std::sort(edges.begin(), edges.end(), me_rank_sort);
  }
  else
  {
    std::sort(edges.begin(), edges.end(), me_score_sort);
  }

  uint64_t stop = 0;

  std::vector<thresh_t> tvec;

  log(LOG_INFO) << "Starting topology assesment\n";
  std::cout.precision(precision);
  for (auto it = steps.rbegin(); it != steps.rend(); it++)
  {
    log(LOG_INFO) << "Current threshold:" << *it << "\n";
    tvec.push_back(check_one(edges, *it, stop, trank));
  }

  get_ranks(tvec, "nn");
  get_ranks(tvec, "ne");
  get_ranks(tvec, "sf");
  get_ranks(tvec, "cc");

  double cut = tvec[0].data["nnr"] + tvec[0].data["sfr"] + tvec[0].data["ccr"];
  thresh_t final = tvec[0];
  for (auto& t : tvec)
  {
    if (t.data["nnr"] + t.data["sfr"] + t.data["ccr"] < cut)
      final = t;
  }

  log(LOG_INFO) << "Suggested threshold: " << final.data["t"] << "("
                << final.data["nnr"] << ',' << final.data["sfr"] << ',' << final.data["ccr"]
                << ")" << '\n';

  rf.close();

  return 0;
}
