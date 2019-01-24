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
#include <tclap/CmdLine.h>
#include <cmath>
#include <boost/numeric/conversion/cast.hpp>

#undef DEBUG

#include <networkit/base/Algorithm.h>
#include <networkit/graph/Graph.h>
#include <networkit/centrality/Centrality.h>
#include <networkit/centrality/EstimateBetweenness.h>
#include <networkit/centrality/ApproxCloseness.h>
#include <networkit/centrality/EigenvectorCentrality.h>
#include <networkit/centrality/KatzCentrality.h>
#include <networkit/centrality/Betweenness.h>
#include <networkit/centrality/PageRank.h>
#include <networkit/centrality/Closeness.h>
#include <networkit/centrality/KPathCentrality.h>
#include <networkit/centrality/SpanningEdgeCentrality.h>
#include <networkit/components/ConnectedComponents.h>

bool in_string(const std::string& query, const std::string& target)
{
  auto p = target.find(query);
  if (p == std::string::npos)
    return false;
  else
    return true;
}

int stats(int argc, char ** argv)
{
  logger log(std::cerr, "stats");

  std::string infile;
  uint32_t tpos;
  bool trank;
  bool exact;
  uint32_t nsamples;
  std::string metrics;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " stats";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Calculate network centrality statistics", ' ', version);

    TCLAP::ValueArg<uint32_t>
    arg_index("i", "index", "Use this index as weights", false,
              0, "float");
    cmd.add(arg_index);

    TCLAP::ValueArg<uint32_t>
    arg_nsamples("n", "nsamples", "Use N samples for approximations", false,
                 0, "float");
    cmd.add(arg_nsamples);

    TCLAP::ValueArg<std::string>
    arg_metrics("m", "metrics", "String describing metrics to calculate"
                " [PR,CLO,BTW,STR,EV,KTZ,SEC,EBC]", false,
                "PR,CLO,BTW,STR,EV,KTZ,SEC,EBC", "string");
    cmd.add(arg_metrics);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true,
               "", "string");
    cmd.add(arg_infile);

    TCLAP::SwitchArg
    arg_r("r", "weight-rank", "Set weight value to rank rather than score",
          cmd, false);

    TCLAP::SwitchArg
    arg_e("e", "exact", "Calculate exact values for centrality", cmd, false);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    tpos = arg_index.getValue();
    trank = arg_r.getValue();
    exact = arg_e.getValue();
    nsamples = arg_nsamples.getValue();
    metrics = arg_metrics.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return EINVAL;
  }

  try
  {
    if (! file_can_read(infile.c_str()))
      throw std::runtime_error("Cannot read input file: " + infile);
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

  if (h.attr.pagerank_calc == 1 && in_string("PR", metrics))
  {
    log(LOG_WARN) << "Found previous PageRank calculation. Deleting...\n";
    h.attr.pagerank_calc = 0;
    h.pagerank.clear();
  }
  if (h.attr.closeness_calc == 1 && in_string("CLO", metrics))
  {
    log(LOG_WARN) << "Found previous closeness calculation. Deleting...\n";
    h.attr.closeness_calc = 0;
    h.closeness.clear();
  }
  if (h.attr.betweenness_calc == 1 && in_string("BTW", metrics))
  {
    log(LOG_WARN) << "Found previous betweenness calculation. Deleting...\n";
    h.attr.betweenness_calc = 0;
    h.betweenness.clear();
  }
  if (h.attr.strength_calc == 1 && in_string("STR", metrics))
  {
    log(LOG_WARN) << "Found previous strength calculation. Deleting...\n";
    h.attr.strength_calc = 0;
    h.strength.clear();
  }
  if (h.attr.eigenvector_calc == 1 && in_string("EV", metrics))
  {
    log(LOG_WARN) << "Found previous eigenvector calculation. Deleting...\n";
    h.attr.eigenvector_calc = 0;
    h.eigenvector.clear();
  }
  if (h.attr.katz_calc == 1 && in_string("KTZ", metrics))
  {
    log(LOG_WARN) << "Found previous katz calculation. Deleting...\n";
    h.attr.katz_calc = 0;
    h.katz.clear();
  }

  if (! exact && nsamples == 0)
  {
    nsamples = boost::numeric_cast<uint64_t>(boost::numeric_cast<double>(h.attr.nodes) * 0.1);
    log(LOG_INFO) << "Using " << nsamples << " as the number of samples\n";
  }

  if (tpos == 0)
    tpos = h.attr.nalgs - 1;
  else
    tpos--;

  log(LOG_INFO) << "Reading graph\n";
  std::vector<SeidrFileEdge> ev =
    read_network(h, rf, -std::numeric_limits<double>::infinity(), tpos, false);

  log(LOG_INFO) << "Creating Networkit graph\n";
  NetworKit::Graph g(h.attr.nodes, true, false);
  for (auto& e : ev)
    g.addEdge(e.index.i, e.index.j, e.scores[tpos].s);

  g.indexEdges();

  log(LOG_INFO) << "Have graph with " << g.numberOfNodes() << " nodes and "
                <<  g.numberOfEdges() << " edges\n";

  rf.close();

  log(LOG_INFO) << "Calculating connected components\n";
  NetworKit::ConnectedComponents ccn(g);
  ccn.run();
  uint64_t comp = ccn.numberOfComponents();

  log(LOG_INFO) << "Have " << comp << " components\n";

  if (in_string("PR", metrics))
  {
    log(LOG_INFO) << "Calculating PageRank\n";
    NetworKit::PageRank prv(g);
    prv.run();
    h.attr.pagerank_calc = 1;
    for (uint32_t i = 0; i < h.attr.nodes; i++)
    {
      h.pagerank.push_back(prv.score(i));
    }
  }

  if (comp == 1 && in_string("CLO", metrics))
  {
    log(LOG_INFO) << "Calculating " << (exact ? "exact" : "approximate")
                  << " Closeness\n";
    NetworKit::ApproxCloseness clv(g, nsamples);
    NetworKit::Closeness clve(g);
    if (exact)
    {
      clve.run();
    }
    else
    {
      clv.run();
    }
    for (uint32_t i = 0; i < h.attr.nodes; i++)
    {
      h.closeness.push_back(exact ? clve.score(i) : clv.score(i));
    }
    h.attr.closeness_calc = 1;
  }
  else
  {
    log(LOG_WARN) << "Graph is disconnected. Closeness is not defined\n";
    for (uint32_t i = 0; i < h.attr.nodes; i++)
    {
      h.closeness.push_back(0);
    }
    h.attr.closeness_calc = 0;
  }

  std::vector<double> ebv;
  if (in_string("BTW", metrics))
  {
    log(LOG_INFO) << "Calculating " << (exact ? "exact" : "approximate")
                  << " Betweenness\n";
    NetworKit::EstimateBetweenness btv(g, nsamples, false, exact);
    NetworKit::Betweenness btve(g, false, true);
    if (exact)
    {
      btve.run();
      if (in_string("EBC", metrics))
      {
        ebv = btve.edgeScores();
      }
    }
    else
    {
      btv.run();
    }
    for (uint32_t i = 0; i < h.attr.nodes; i++)
    {
      h.betweenness.push_back(exact ? btve.score(i) : btv.score(i));
    }
    h.attr.betweenness_calc = 1;
  }

  if (in_string("EV", metrics))
  {
    log(LOG_INFO) << "Calculating Eigenvector centrality\n";
    NetworKit::EigenvectorCentrality evv(g);
    evv.run();
    h.attr.eigenvector_calc = 1;
    for (uint32_t i = 0; i < h.attr.nodes; i++)
    {
      h.eigenvector.push_back(evv.score(i));
    }
  }

  if (in_string("KTZ", metrics))
  {
    log(LOG_INFO) << "Calculating Katz centrality\n";
    NetworKit::KatzCentrality kav(g);
    kav.run();
    h.attr.katz_calc = 1;
    for (uint32_t i = 0; i < h.attr.nodes; i++)
    {
      h.katz.push_back(kav.score(i));
    }
  }

  if (in_string("STR", metrics))
  {
    log(LOG_INFO) << "Calculating Node Strength\n";
    h.attr.strength_calc = 1;
    for (uint64_t i = 0; i < h.attr.nodes; i++)
    {
      h.strength.push_back(g.weightedDegree(i));
    }
  }

  std::vector<double> spv;
  if (in_string("SEC", metrics))
  {
    log(LOG_INFO) << "Calculating Spanning Edge centrality\n";

    NetworKit::SpanningEdgeCentrality spc(g);
    if (exact)
      spc.run();
    else
      spc.runParallelApproximation();
    spv = spc.scores(true);
  }

  log(LOG_INFO) << "Writing results\n";


  bool sec_found = false;
  bool ebc_found = false;
  size_t sec_i = 0;
  size_t ebc_i = 0;

  if (in_string("SEC", metrics))
  {
    auto sec_ptr = std::find(h.supp.begin(), h.supp.end(), "SEC");
    if (sec_ptr == h.supp.end())
    {
      h.attr.nsupp += 1;
      h.attr.nsupp_flt += 1;
      h.supp.push_back("SEC");
    }
    else
    {
      sec_found = true;
      sec_i = std::distance(h.supp.begin(), sec_ptr);
      sec_i -= h.attr.nsupp_str;
      sec_i -= h.attr.nsupp_int;
    }
  }

  if (in_string("EBC", metrics))
  {
    auto ebc_ptr = std::find(h.supp.begin(), h.supp.end(), "EBC");
    if (ebc_ptr == h.supp.end())
    {
      h.attr.nsupp += 1;
      h.attr.nsupp_flt += 1;
      h.supp.push_back("EBC");
    }
    else
    {
      ebc_found = true;
      ebc_i = std::distance(h.supp.begin(), ebc_ptr);
      ebc_i -= h.attr.nsupp_str;
      ebc_i -= h.attr.nsupp_int;
    }
  }

  std::string tempfile = infile + ".tmp";
  SeidrFile tf(tempfile.c_str());
  tf.open("w");

  h.serialize(tf);

  for (uint64_t i = 0; i < ev.size(); i++)
  {
    SeidrFileEdge x = ev[i];
    if (in_string("SEC", metrics))
    {
      if (sec_found)
        x.supp_flt[sec_i] = spv[i];
      else
        x.supp_flt.push_back(spv[i]);
    }
    if (in_string("EBC", metrics))
    {
      if (ebc_found)
        x.supp_flt[ebc_i] = ebv[i];
      else
        x.supp_flt.push_back(ebv[i]);
    }
    x.serialize(tf, h);
  }

  tf.close();

  rename(tempfile, infile);

  return 0;
}
