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

// Seir
#include <common.h>
#include <fs.h>
#include <BSlogger.h>
#include <Serialize.h>
#include <backbone.h>
// External
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

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

  seidr_backbone_param_t param;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " backbone";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella("Determine noisy network backbone scores. "
                                   "Optionally filter on these scores");

  po::options_description opt("Common Options");
  opt.add_options()
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("help,h", "Show this help message")
  ("out-file,o", po::value<std::string>(&param.outfile)->default_value("auto"),
   "Output file name ['-' for stdout]")
  ("tempdir,T",
   po::value<std::string>(&param.tempdir)->default_value("", "auto"),
   "Directory to store temporary data");

  po::options_description bopt("Backbone Options");
  bopt.add_options()
  ("index,i", po::value<uint32_t>(&param.tpos)->default_value(0, "last index"),
   "Score index to use")
  ("filter,F", po::value<double>(&param.filter)->default_value(-1, "no filter"),
   "Subset network to edges with at least this SD. "
   "1.28, 1.64, and 2.32 correspond to ~P0.01, 0.05 and 0.1.");

  po::options_description req("Required [can be positional]");
  req.add_options()
  ("in-file", po::value<std::string>(&param.infile)->required(),
   "Input SeidrFile");

  umbrella.add(req).add(bopt).add(opt);

  po::positional_options_description p;
  p.add("in-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, args).
            options(umbrella).positional(p).run(), vm);

  if (vm.count("help") || argc == 1)
  {
    std::cerr << umbrella << '\n';
    return EINVAL;
  }

  try
  {
    po::notify(vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  if (param.outfile == "auto")
  {
    param.outfile = to_absolute(param.infile);
    param.outfile = replace_ext(param.outfile, "bb.sf");
  }

  param.tempfile = tempfile(param.tempdir);

  try
  {
    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);
    assert_dir_is_writeable(dirname(param.outfile));
    if (! param.force)
    {
      assert_no_overwrite(param.outfile);
    }
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  SeidrFile rf(param.infile.c_str());
  rf.open("r");

  SeidrFileHeader h;
  h.unserialize(rf);

  make_tpos(param.tpos, h);

  uint32_t cnt = 0;
  std::map<std::string, uint32_t> node_map;
  for (std::string& node : h.nodes)
    node_map[node] = cnt++;
  rf.close();

  // Setup phase. Calculate global sum of edge weight and per-node sum of
  // edge weights
  double global_strength = 0;
  std::vector<double> node_strengths;
  node_strengths.resize(h.nodes.size(), 0);
  bool nc_calculcated = h.have_supp("NC_Score");
  if (! nc_calculcated)
  {
    log(LOG_INFO) << "First pass. Getting global node strengths\n";
    // Reset reading position to beginning of file
    rf.open("r");
    rf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
      if (std::isfinite(e.scores[param.tpos].s))
      {
        global_strength += e.scores[param.tpos].s;
        node_strengths[e.index.i] += e.scores[param.tpos].s;
        node_strengths[e.index.j] += e.scores[param.tpos].s;
      }
    });
    rf.close();
  }

  // Add scores for upper triangle
  global_strength *= 2;
  // for(auto& s : node_strengths)
  //   s *= 2;

  // Reset reading position to beginning of file
  rf.open("r");

  // Set up output
  log(LOG_INFO) << "Using temp file: " << param.tempfile << '\n';
  SeidrFile out(param.tempfile.c_str());
  out.open("w");
  if (! nc_calculcated)
  {
    log(LOG_INFO) << "Calculating NC Score and NC SDev...\n";
    h.attr.nsupp += 2;
    h.attr.nsupp_flt += 2;
    h.supp.push_back("NC_Score");
    h.supp.push_back("NC_SDev");
  }
  h.serialize(out);

  uint16_t nc_score_ind = 0;
  uint16_t nc_sdev_ind = 0;

  if (nc_calculcated)
  {
    log(LOG_INFO) << "Found already calculated NC Score and NC SDev."
                  " Using those\n";
    nc_score_ind = h.get_supp_ind("NC_Score");
    nc_sdev_ind = h.get_supp_ind("NC_SDev");
  }
  else
  {
    nc_score_ind = h.attr.nsupp - h.attr.nsupp_str - h.attr.nsupp_int - 2;
    nc_sdev_ind = nc_score_ind + 1;
  }

  uint64_t ectr = 0;

  // Main phase, calculate P-value based on binomial cdf
  rf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & hs) {
    if (! nc_calculcated) // Need to add the backbone data
    {
      if (std::isfinite(e.scores[param.tpos].s))
      {
        double& ni = node_strengths[e.index.i];
        double& nj = node_strengths[e.index.j];
        double nij = e.scores[param.tpos].s;
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
    }
    if (param.filter > 0)
    {
      if (e.supp_flt[nc_score_ind] - (e.supp_flt[nc_sdev_ind] * param.filter) > 0)
      {
        e.serialize(out, h);
        ectr++;
      }
      else if (h.attr.dense)
      {
        SeidrFileEdge ex;
        EDGE_SET_MISSING(ex.attr.flag);
        ex.serialize(out, h);
      }
    }
    else
    {
      e.serialize(out, h);
      ectr++;
    }

  });

  rf.close();
  out.close();

  log(LOG_INFO) << "Done. Edges deleted: " << h.attr.edges - ectr << '\n';
  log(LOG_INFO) << "Edges remaining: " << ectr << '\n';

  h.attr.edges = ectr;
  double nn = h.attr.nodes;
  double ne = ectr;
  if (ne > ((nn * (nn - 1) / 2) * 0.66))
    h.attr.dense = 1;
  else
    h.attr.dense = 0;

  SeidrFile rfx(param.outfile.c_str());
  rfx.open("w");
  h.serialize(rfx);

  out.open("r");
  SeidrFileHeader hx;
  hx.unserialize(out);

  log(LOG_INFO) << "Saving new file...\n";
  if (hx.attr.dense)
  {
    for (uint64_t i = 0; i < hx.attr.nodes; i++)
    {
      for (uint64_t j = 0; j < i; j++)
      {
        SeidrFileEdge e;
        e.unserialize(out, hx);
        e.index.i = i;
        e.index.j = j;
        if (h.attr.dense)
        {
          e.serialize(rfx, h);
        }
        else if (EDGE_EXISTS(e.attr.flag))
        {
          e.serialize(rfx, h);
        }
      }
    }
  }
  else
  {
    for (uint64_t i = 0; i < ectr; i++)
    {
      SeidrFileEdge e;
      e.unserialize(out, hx);
      e.serialize(rfx, h);
    }
  }

  out.close();
  rfx.close();

  remove(param.tempfile);

  if (param.filter > 0)
  {
    log(LOG_WARN) << "Unless you need to keep node indices consistent, you "
                  << "probably want to run 'seidr reheader' "
                  << "on the filtered file\n";
  }

  return 0;
}
