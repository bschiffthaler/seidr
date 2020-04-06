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
#include <common.h>
#include <tom_fun.h>
#include <fs.h>

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

#include <armadillo>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char ** argv)
{
  logger log(std::cerr, "TOMSimilarity");

  seidr_tom_param_t param;

  po::options_description umbrella("Topological overlap matrix as implemented "
                                   "in WGCNA");
  po::variables_map vm;
  try
  {
    po::options_description opt("Common Options");
    opt.add_options()
    ("help,h", "Show this help message")
    ("targets,t", po::value<std::string>(&param.targets_file),
     "File containing gene names"
     " of genes of interest. The network will only be"
     " calculated using these as the sources of potential connections.")
    ("outfile,o",
     po::value<std::string>(&param.outfile)->default_value("tom_scores.tsv"),
     "Output file path")
    ("verbosity,v",
     po::value<unsigned>(&param.verbosity)->default_value(3),
     "Verbosity level (lower is less verbose)")
    ("force,f", po::bool_switch(&param.force)->default_value(false),
     "Force overwrite if output already exists");

    po::options_description algopt("Correlation specific Options");
    algopt.add_options()
    ("method,m",
     po::value<std::string>(&param.method)->default_value("pearson"),
     "Correlation method: [bicor, pearson]")
    ("absolute,a", po::bool_switch(&param.abs)->default_value(false),
     "Report absolute values")
    ("scale,s", po::bool_switch(&param.do_scale)->default_value(false),
     "Transform data to z-scores")
    ("sft,b", po::value<uint64_t>(&param.sft)->default_value(0),
     "Soft threshold to apply (0 for autodetection)")
    ("max-power,M", po::value<uint64_t>(&param.max_power)->default_value(30),
     "Maximum power to check for SFT test in auto detection mode")
    ("sft-cutoff,S", po::value<double>(&param.sft_cutoff)->default_value(0.8, "0.8"),
     "Soft threshold R^2 cutoff in autodetection mode")
    ("tom-type,T",
     po::value<std::string>(&param.tom_type)->default_value("signed"),
     "TOM type: unsigned, signed, signed-hybrid")
    ;

    po::options_description req("Required Options");
    req.add_options()
    ("infile,i", po::value<std::string>(&param.infile)->required(),
     "The expression table (without headers)")
    ("genes,g", po::value<std::string>(&param.gene_file)->required(),
     "File containing gene names");

    umbrella.add(req).add(algopt).add(opt);

    po::positional_options_description p;


    po::store(po::command_line_parser(argc, argv).
              options(umbrella).run(), vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << "Argument exception: " << e.what() << '\n';
    return 22;
  }


  if (vm.count("help") > 0 || argc == 1)
  {
    std::cerr << umbrella << '\n';
    return 1;
  }

  try
  {
    po::notify(vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << "Argument exception: " << e.what() << '\n';
    return 22;
  }

  log.set_log_level(param.verbosity);

  try
  {
    param.outfile = to_absolute(param.outfile);
    param.infile = to_absolute(param.infile);
    param.gene_file = to_absolute(param.gene_file);

    assert_exists(dirname(param.outfile));
    assert_exists(param.infile);
    assert_is_regular_file(param.infile);
    assert_exists(param.gene_file);
    assert_can_read(param.gene_file);
    assert_can_read(param.infile);

    if (vm.count("targets") > 0)
    {
      param.targets_file = to_absolute(param.targets_file);
      assert_exists(param.targets_file);
    }

    if (vm["outfile"].defaulted())
    {
      param.outfile = param.method + "_tom_scores.tsv";
    }

    if (! param.force)
    {
      assert_no_overwrite(param.outfile);
    }
    assert_arg_constraint<std::string>({"pearson", "bicor"}, param.method);
    assert_arg_constraint<std::string>({"unsigned", "signed", "signed-hybrid"},
                                       param.tom_type);
    assert_in_range<double>(param.sft_cutoff, 0, 1, "sft-cutoff");
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  arma::mat gm;

  try
  {
    gm.load(param.infile, arma::raw_ascii);
    verify_matrix(gm);
    if (param.do_scale)
    {
      log(LOG_INFO) << "Transforming matrix to z-score\n";
      scale(gm);
    }

    if (param.method == "pearson")
    {
      gm = arma::cor(gm);
    }
    else if (param.method == "bicor")
    {
      gm = bicor(gm);
    }

    uint64_t p = param.sft;
    tom_cor_t tt;
    if (param.tom_type == "unsigned")
    {
      tt = UNSIGNED;
    }
    else if (param.tom_type == "signed")
    {
      tt = SIGNED;
    }
    else if (param.tom_type == "signed-hybrid")
    {
      tt = SIGNED_HYBRID;
    }
    if (p == 0)
    {
      for (uint64_t i = 2; i <= param.max_power; i += 2)
      {
        log(LOG_INFO) << "Checking scale free fit. SFT = " << i << '\n';
        double sft = scale_free_fit(gm, i, 10, tt);
        p = i;
        log(LOG_INFO) << "SFT (" << i << ") = " << sft << '\n';
        if (sft > 0.8)
        {
          log(LOG_INFO) << "Power " << i << " satisfied fit.\n";
          break;
        }
      }
    }

    gm = tom_similarity(gm, p, tt);

    if (vm.count("targets") == 0)
    {
      write_lm(gm, param.outfile, param.abs);
    }
    else
    {
      std::vector<std::string> genes = read_genes(param.gene_file);
      std::vector<std::string> targets = read_genes(param.targets_file);
      std::unordered_map<std::string, uint64_t> gene_map;
      uint64_t ctr = 0;
      for (auto& g : genes)
      {
        gene_map[g] = ctr++;
      }
      std::ofstream ofs(param.outfile, std::ios::out);
      for (auto& t : targets)
      {
        uint64_t i;
        try
        {
          i = gene_map.at(t);
        }
        catch (std::exception& e)
        {
          log(LOG_ERR) << e.what() << '\n';
          log(LOG_ERR) << "Target gene " << t << "is not in the expression matrix\n";
        }
        for (uint64_t j = 0; j < gm.n_cols; j++)
        {
          if (i == j)
          {
            continue;
          }
          ofs << genes[i] << '\t' << genes[j] << '\t' <<
              (param.abs ? fabs(gm(i, j)) : gm(i, j)) << '\n';
        }
      }
    }
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }
  return 0;
}
