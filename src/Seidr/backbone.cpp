// Seir
#include <common.h>
#include <fs.h>
#include <BSlogger.h>
#include <Serialize.h>
#include <backbone.h>
#include <tclap/CmdLine.h>
// External
#include <cmath>
#include <map>
#include <vector>
#include <string>

// Binomial cumulative distribution function adjusted for seidr where
// k always [0, 1]
double binom_cdf(const double& k, const double& n, const double& p)
{
  double fin = 0;
  // Case ki < 1
  fin += std::pow(1 - p, n);
  // Special case for ki == 1 (the best score)
  if ( almost_equal(k, 1) )
  {
    fin += (n * p * std::pow(1 - p, n - 1));
  }
  return fin;
}

double prior_prob(const double& ni, const double& nj, const double& n)
{
  return ((ni * nj) / n) * (1 / n);
}

int backbone(int argc, char * argv[])
{

  logger log(std::cerr, "backbone");

  std::string infile;
  uint32_t tpos;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " neighbours";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Calculate P values to determine network backbone.", ' ',
                       version);

    TCLAP::ValueArg<uint32_t>
    arg_tpos("i", "index", "Index of score to use", false,
             0, "last column");
    cmd.add(arg_tpos);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true, "", "file");
    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    tpos = arg_tpos.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  try
  {
    infile = to_canonical(infile);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    log(LOG_ERR) << "Make sure that the file " << infile << " exists.\n";
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

  uint32_t cnt = 0;
  std::map<std::string, uint32_t> node_map;
  for (std::string& node : h.nodes)
    node_map[node] = cnt++;
  rf.close();

  // Reset reading position to beginning of file
  rf.open("r");

  // Setup phase. Calculate global sum of edge weight and per-node sum of
  // edge weights
  double global_strength = 0;
  std::vector<double> node_strengths;
  node_strengths.resize(h.nodes.size(), 0);
  rf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    if (std::isfinite(e.scores[tpos].s))
    {
      global_strength += e.scores[tpos].s;
      node_strengths[e.index.i] += e.scores[tpos].s;
      node_strengths[e.index.j] += e.scores[tpos].s;
    }
  });
  rf.close();

  // Add scores for upper triangle
  global_strength *= 2;
  // for(auto& s : node_strengths)
  //   s *= 2;

  // Reset reading position to beginning of file
  rf.open("r");

  // Set up output
  std::string tmp = infile + ".tmp";
  SeidrFile out(tmp.c_str());
  out.open("w");
  h.attr.nsupp += 2;
  h.attr.nsupp_flt += 2;
  h.supp.push_back("NC_Score");
  h.supp.push_back("NC_SDev");
  h.serialize(out);

  // Main phase, calculate P-value based on binomial cdf
  rf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & hs) {
    if (std::isfinite(e.scores[tpos].s))
    {
      double& ni = node_strengths[e.index.i];
      double& nj = node_strengths[e.index.j];
      double& nij = e.scores[tpos].s;
      // Essentially a verbatim copy of Coscia's & Neffe's 
      // (https://ieeexplore.ieee.org/document/7929996/) python source. If you
      // need this in actually readable form, I suggest chapter IV of the 
      // aforementioned paper 
      double prior = prior_prob(ni, nj, global_strength);
      double kappa = global_strength / (ni * nj);
      double score = ((kappa * nij) - 1) / ((kappa * nij + 1));
      double var_prior_prob = (1 / (global_strength * global_strength)) *
                              (ni * nj *
                               (global_strength - ni) *
                               (global_strength - nj)) /
                              ((global_strength * global_strength) *
                               ((global_strength - 1)));
      double alpha_prior = (((prior * prior) /
                             var_prior_prob) * (1 - prior)) - prior;
      double beta_prior = (prior / var_prior_prob) * (1 - (prior * prior)) -
                          (1 - prior);
      double alpha_post = alpha_prior + nij;
      double beta_post = global_strength - nij + beta_prior;
      double expected_pij = alpha_post / (alpha_post + beta_post);
      double variance_nij = expected_pij * (1 - expected_pij) * global_strength;
      double d = (1 / (ni * nj)) -
                 (global_strength * ((ni + nj) / (std::pow((ni * nj), 2))));
      double variance_cij = variance_nij * std::pow(((2 * (kappa + (nij * d))) /
                            ( std::pow((kappa * nij) + 1, 2))), 2);
      double sdev_cij = std::pow(variance_cij, 0.5);

      e.supp_flt.push_back(score);
      e.supp_flt.push_back(sdev_cij);
    }
    else
    {
      e.supp_flt.push_back(NAN);
      e.supp_flt.push_back(NAN);
    }
    e.serialize(out, h);
  });

  out.close();

  return 0;

}
