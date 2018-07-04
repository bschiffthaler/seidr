#include <common.h>
#include <stats.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <cmath>
#include <boost/numeric/conversion/cast.hpp>

#undef DEBUG

#include <networkit/base/Algorithm.h>
#include <networkit/graph/Graph.h>
#include <networkit/distance/APSP.h>
#include <networkit/components/ConnectedComponents.h>
#include <networkit/global/ClusteringCoefficient.h>

std::pair<double, double> SFT(const NetworKit::Graph& g, uint64_t n_bins = 200)
{
  uint64_t nn = g.numberOfNodes();
  double nnd = boost::numeric_cast<double>(nn);

  arma::vec deg_vec(nn, arma::fill::zeros);

  for(uint64_t i = 0; i < nn; i++)
  {
    deg_vec(i) = g.degree(i);
  }

  double deg_sum = arma::as_scalar(arma::accu(deg_vec));
  double avg_deg = deg_sum / nnd;

  arma::vec centers = arma::linspace<arma::vec>
  (arma::min(deg_vec), arma::max(deg_vec), n_bins);

  arma::vec hist = arma::conv_to<arma::vec>::from(arma::hist(deg_vec, centers));

  centers = arma::log10(centers);
  hist = arma::log10(hist + 1); //Add 1 to avoid log10(0)

  // For the most part R^2 of a linear model is equal to the correlation squared
  double rsq = arma::as_scalar(arma::pow(arma::cor(centers, hist), 2));

  return std::pair<double, double>(avg_deg, rsq);
}

double avg_w_deg(const NetworKit::Graph& g)
{

  uint64_t nn = g.numberOfNodes();
  double nnd = boost::numeric_cast<double>(nn);

  arma::vec deg_vec(nn, arma::fill::zeros);

  for(uint64_t i = 0; i < nn; i++)
  {
    deg_vec(i) = g.weightedDegree(i);
  }

  double deg_sum = arma::as_scalar(arma::accu(deg_vec));
  double avg_deg = deg_sum / nnd;

  return avg_deg;
}

std::pair<double, double> 
diameter(const std::vector<std::vector<double>>& sps)
{
  double max = -std::numeric_limits<double>::infinity();
  double nn = boost::numeric_cast<double>(sps.size());
  double sum = 0;

  for(const auto& node : sps)
  {
    for(const auto& pathlength : node)
    {
      sum += pathlength;
      if (pathlength > max) max = pathlength;
    }
  }

  return std::pair<double, double>(max, sum / nn);
}

int graphstats(int argc, char ** argv)
{
  logger log(std::cerr, "graphstats");
  
  std::string infile;
  uint32_t tpos;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " graphstats";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Calculate overall network statistics", ' ', version);

    TCLAP::ValueArg<uint32_t>
    arg_index("i", "index", "Use this index as weights", false,
              0, "float");
    cmd.add(arg_index);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true,
               "", "string");
    cmd.add(arg_infile);

    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    tpos = arg_index.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return EINVAL;
  }

  SeidrFile rf(infile.c_str());
  rf.open("r");

  SeidrFileHeader h;
  h.unserialize(rf);

  log(LOG_INFO) << "Reading graph\n";
  std::vector<SeidrFileEdge> ev =
    read_network(h, rf, -std::numeric_limits<double>::infinity(), tpos, false);

  log(LOG_INFO) << "Creating Networkit graph\n";
  NetworKit::Graph g(h.attr.nodes, true, false);
  for (auto& e : ev)
    g.addEdge(e.index.i, e.index.j, e.scores[tpos].s);

  std::cout << "Number of Nodes:\t" << g.numberOfNodes() << '\n';
  std::cout << "Number of Edges:\t" << g.numberOfEdges() << '\n';

  log(LOG_INFO) << "Calculating connected components\n";
  NetworKit::ConnectedComponents ccn(g);
  ccn.run();
  uint64_t comp = ccn.numberOfComponents();
  std::cout << "Number of Connected Components:\t" << comp << '\n';

  log(LOG_INFO) << "Calculating clustering coefficient\n";
  NetworKit::ClusteringCoefficient ccoef;
  double gccf = ccoef.exactGlobal(g);
  std::cout << "Global clustering coefficient:\t" << gccf << '\n';

  log(LOG_INFO) << "Calculating scale-freeness\n";
  auto sft = SFT(g);
  std::cout << "Scale free fit:\t" << sft.second << '\n';
  std::cout << "Average degree:\t" << sft.first << '\n';
  std::cout << "Average weighted degree:\t" << avg_w_deg(g) << '\n';

  NetworKit::APSP apsp(g);
  apsp.run();
  auto all_dist = apsp.getDistances();
  auto dia = diameter(all_dist);
  std::cout << "Network diameter:\t" << dia.first << '\n';
  std::cout << "Average path length:\t" << dia.second << '\n';

  return 0;
}