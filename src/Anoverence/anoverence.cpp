// Seidr
#include <common.h>
#include <fs.h>
#include <anova-fun.h>
#include <BSlogger.h>
// External
#include <string>
#include <vector>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#define ANOVA_FULL 0
#define ANOVA_PARTIAL 1

namespace fs = boost::filesystem;

int main(int argc, char ** argv) {

  std::string infile;
  float weight;
  std::string feature_file;
  size_t mode = ANOVA_FULL;
  std::string outfile;
  bool force = false;
  std::string targets_file;
  std::string gene_file;

  LOG_INIT_CERR();

  // Define program options
  try
  {
    TCLAP::CmdLine cmd("ANOVERENCE implementation for seidr", ' ', version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<double>
    weight_arg("w", "weight", "Weight for knockout genes",
               true, 1.0, "1.0");
    cmd.add(weight_arg);

    TCLAP::ValueArg<std::string>
    feature_arg("e", "features", "File containing experiment metadata.", true,
                "", "");
    cmd.add(feature_arg);

    TCLAP::ValueArg<std::string>
    outfile_arg("o", "outfile", "Output file path", false, "edgelist.tsv",
                "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::ValueArg<std::string>
    targets_arg("t", "targets", "File containing target gene ids", false,
                "", "");
    cmd.add(targets_arg);

    TCLAP::ValueArg<std::string>
    genes_arg("g", "genes", "File containing gene ids", false, "", "");
    cmd.add(genes_arg);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.parse(argc, argv);

    infile = infile_arg.getValue();
    outfile = outfile_arg.getValue();
    feature_file = feature_arg.getValue();
    weight = weight_arg.getValue();
    targets_file = targets_arg.getValue();
    gene_file = genes_arg.getValue();
    force = switch_force.getValue();

    if (targets_file != "")
      mode = ANOVA_PARTIAL;

  }
  catch (TCLAP::ArgException& e)  // catch any exceptions
  {
    log(LOG_ERR) << e.error() << " for arg " << e.argId() << '\n';
    return EINVAL;
  }

  // Check all kinds of FS problems that may arise
  try
  {
    outfile = to_absolute(outfile);
    infile = to_absolute(infile);

    if (! file_exists(dirname(outfile)) )
      throw std::runtime_error("Directory does not exist: " +
                               dirname(outfile));

    if (! file_exists(infile) )
      throw std::runtime_error("File does not exist: " + infile);

    if (! regular_file(infile) )
      throw std::runtime_error("Not a regular file: " + infile);

    if (! file_can_read(infile) )
      throw std::runtime_error("Cannot read: " + infile);

    if (weight < 0 || almost_equal(weight, 0))
      throw std::runtime_error("Weight can't be less than or equal 0");


    if (mode == ANOVA_PARTIAL)
    {
      targets_file = to_absolute(targets_file);
      if (! file_exists(targets_file) )
        throw std::runtime_error("File does not exist: " + targets_file);

      if (! file_can_read(targets_file))
        throw std::runtime_error("Cannot read: " + targets_file);

      gene_file = to_absolute(gene_file);
      if (! file_exists(gene_file) )
        throw std::runtime_error("File does not exist: " + gene_file);

      if (! file_can_read(gene_file))
        throw std::runtime_error("Cannot read: " + gene_file);

    }

    if (! force && file_exists(outfile))
      throw std::runtime_error("File exists: " + outfile);

  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return EINVAL;
  }

  try
  {
    // Get input files
    readExpressionData(infile);
    readFeatures(feature_file);
    readGeneFile(gene_file);

    std::vector<std::string> targets;
    if (mode == ANOVA_PARTIAL)
      targets = read_genes(targets_file, '\t', '\n');


    // Common functions that are always called
    mapChipVector();
    calculateMean();
    selectPairs();


    switch (mode) {
    case ANOVA_FULL:
      analizeGenePairs(outfile, weight);
      break;

    case ANOVA_PARTIAL:
      analizePartialPairs(targets, outfile, weight);
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
