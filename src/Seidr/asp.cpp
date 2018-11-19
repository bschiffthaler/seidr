#include <common.h>
#include <asp.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>

#undef DEBUG

#include <networkit/graph/Graph.h>
#include <networkit/distance/APSP.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <tclap/CmdLine.h>
#include <cmath>
#include <boost/numeric/conversion/cast.hpp>



int asp(int argc, char ** argv)
{
  logger log(std::cerr, "asp");

  uint32_t tpos;
  std::string infile;
  bool trank;
  uint16_t precision;
  bool invert;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " asp";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Compute all shortest paths in the network", ' ',
                       version);

    TCLAP::ValueArg<uint32_t>
    arg_index("i", "index", "Use this index", false,
              0, "last column");
    cmd.add(arg_index);

    TCLAP::ValueArg<uint16_t>
    arg_precision("p", "precision", "Number of digits after comma", false,
                  8, "8");
    cmd.add(arg_precision);

    TCLAP::SwitchArg
    arg_r("r", "use-rank",
          "Use rank rather than score", cmd, false);

    TCLAP::SwitchArg
    arg_invert("I", "invert",
               "Invert rank/score", cmd, false);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true,
               "", "");

    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    tpos = arg_index.getValue();
    trank = arg_r.getValue();
    precision = arg_precision.getValue();
    invert = arg_invert.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return EINVAL;
  }

  try
  {
    file_can_read(infile.c_str());
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

  log(LOG_INFO) << "Reading network\n";
  std::vector<MiniEdge> edges =
    read_network_minimal(h, rf, -std::numeric_limits<double>::infinity(),
                         tpos, trank);
  log(LOG_INFO) << "Read " << edges.size() << " edges\n";

  NetworKit::Graph g(h.attr.nodes, true, false);
  for (uint64_t i = 0; i < edges.size(); i++)
  {
    double weight = edges[i].s;
    if (almost_equal(weight, 0))
      weight = 1e-8;
    weight = invert ? 1.0 / weight : weight;
    g.addEdge(edges[i].i, edges[i].j, weight);
  }

  NetworKit::APSP apsp(g);
  apsp.run();

  for (uint64_t i = 1; i < h.attr.nodes; i++)
  {
    for (uint64_t j = 0; j < i; j++)
    {
      std::cout << apsp.getDistance(i, j) 
      << ((j == i - 1) ? '\n' : '\t');
    }
  }

  rf.close();

  return 0;
}
