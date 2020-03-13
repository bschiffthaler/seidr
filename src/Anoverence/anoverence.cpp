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

// Seidr
#include <common.h>
#include <fs.h>
#include <anova-fun.h>
#include <BSlogger.hpp>
// External
#include <string>
#include <vector>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char ** argv) {

  LOG_INIT_CERR();

  seidr_anova_param_t param;

  po::options_description umbrella("ANOVERENCE implementation for Seidr");

  po::options_description opt("Common Options");
  opt.add_options()
  ("help,h", "Show this help message")
  ("targets,t", po::value<std::string>(&param.targets_file),
   "File containing gene names"
   " of genes of interest. The network will only be"
   " calculated using these as the sources of potential connections.")
  ("outfile,o",
   po::value<std::string>(&param.outfile)->default_value("anova_scores.tsv"),
   "Output file path")
  ("verbosity,v",
   po::value<unsigned>(&param.verbosity)->default_value(3),
   "Verbosity level (lower is less verbose)")
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite if output already exists");

  po::options_description algopt("ANOVERENCE specific Options");
  algopt.add_options()
  ("weight,w",
   po::value<float>(&param.weight)->default_value(1.0, "1.0"),
   "Weight for knockout genes");

  po::options_description req("Required Options");
  req.add_options()
  ("infile,i", po::value<std::string>(&param.infile)->required(),
   "The expression table (without headers)")
  ("genes,g", po::value<std::string>(&param.gene_file)->required(),
   "File containing gene names")
  ("features,e", po::value<std::string>(&param.feature_file)->required(),
   "File containing experiment metadata");

  umbrella.add(req).add(algopt).add(opt);

  po::positional_options_description p;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(umbrella).run(), vm);

  if (vm.count("help") || argc == 1)
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
  }

  log.set_log_level(param.verbosity);

  try
  {
    param.outfile = to_absolute(param.outfile);
    param.infile = to_absolute(param.infile);
    param.gene_file = to_absolute(param.gene_file);

    assert_exists(dirname(param.outfile));
    assert_exists(param.infile);
    assert_exists(param.feature_file);
    assert_is_regular_file(param.infile);
    assert_exists(param.gene_file);
    assert_can_read(param.gene_file);
    assert_can_read(param.infile);
    assert_can_read(param.feature_file);

    if (vm.count("targets"))
    {
      param.targets_file = to_absolute(param.targets_file);
      assert_exists(param.targets_file);
    }

    if (! param.force)
      assert_no_overwrite(param.outfile);

  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  try
  {
    // Get input files
    int geneNamesSize = readGeneFile(param.gene_file);
    int geneNumber = readExpressionData(param.infile);

    if (geneNamesSize != geneNumber) {
      throw std::runtime_error("Data genes and gene names total are not equal.\n");
    }
    readFeatures(param.feature_file);

    if (vm.count("targets"))
      param.mode = ANOVA_PARTIAL;


    std::vector<std::string> targets;
    if (param.mode == ANOVA_PARTIAL)
      targets = read_genes(param.targets_file, '\t', '\n');


    // Common functions that are always called
    mapChipVector();
    calculateMean();
    selectPairs();


    switch (param.mode) {
    case ANOVA_FULL:
      analizeGenePairs(param.outfile, param.weight);
      break;

    case ANOVA_PARTIAL:
      analizePartialPairs(targets, param.outfile, param.weight);
      break;

    default:
      return 1;
    }
  }
  catch (const std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }
  catch (const std::exception& e)
  { //TODO treat other exceptions differently
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }
  return 0;

}
