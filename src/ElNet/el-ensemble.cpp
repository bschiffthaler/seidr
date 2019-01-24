// Seidr
#include <common.h>
#include <fs.h>
#include <mpims.h>
#include <elnet-fun.h>
// External
#include <mpi.h>
#include <svm.h>
#include <cerrno>
#include <string>
#include <vector>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

#define EL_FULL 0
#define EL_PARTIAL 1

namespace fs = boost::filesystem;

int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  seidr_mpi_logger log;

  arma::mat gene_matrix;
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  size_t bs;
  size_t mode = EL_FULL;
  std::string outfile;
  std::vector<std::string> genes;
  std::vector<std::string> targets;
  arma::uword min_sample_size;
  arma::uword max_sample_size;
  arma::uword predictor_sample_size_min;
  arma::uword predictor_sample_size_max;
  arma::uword ensemble_size;
  double alpha;
  double flmin;
  arma::uword nlam;
  unsigned verbosity;

  // Define program options
  try
  {
    TCLAP::CmdLine cmd("Elastic net ensemble implementation for Seidr", ' ',
                       version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<std::string>
    genefile_arg("g", "genes", "File containing gene names", true, "",
                 "string");
    cmd.add(genefile_arg);

    TCLAP::ValueArg<std::string>
    targets_arg("t", "targets", "File containing gene names"
                " of genes of interest. The network will only be"
                " calculated using these as the sources of potential connections.",
                false, "", "string");
    cmd.add(targets_arg);

    TCLAP::ValueArg<unsigned long>
    batchsize_arg("b", "batch-size", "Number of genes in MPI batch", false, 20,
                  "20");
    cmd.add(batchsize_arg);

    TCLAP::ValueArg<std::string> outfile_arg("o", "outfile", "Output file path",
        false, "edgelist.tsv", "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::ValueArg<seidr_uword_t> ensemble_arg("e", "ensemble",
        "The ensemble size", false, 1000, "1000");
    cmd.add(ensemble_arg);

    TCLAP::ValueArg<seidr_uword_t> min_pred_arg("p", "min-predictor-size",
        "The minimum number of predictors to be sampled.", false,
        0, "1/5th of genes");
    cmd.add(min_pred_arg);

    TCLAP::ValueArg<seidr_uword_t> max_pred_arg("P", "max-predictor-size",
        "The maximum number of preditors to be sampled", false,
        0, "4/5th of genes");
    cmd.add(max_pred_arg);

    TCLAP::ValueArg<seidr_uword_t> min_exp_arg("x", "min-experiment-size",
        "The minimum number of experiments to be sampled.", false,
        0, "1/5th of experiments");
    cmd.add(min_exp_arg);

    TCLAP::ValueArg<seidr_uword_t> max_exp_arg("X", "max-experiment-size",
        "The maximum number of experiments to be sampled", false,
        0, "4/5th of experiments");
    cmd.add(max_exp_arg);

    TCLAP::ValueArg<seidr_uword_t>
    nlam_arg("n", "nlambda", "The maximum number of lambda values", false,
             10, "10");
    cmd.add(nlam_arg);

    TCLAP::ValueArg<double>
    alpha_arg("a", "alpha", "The elastic net mixing value alpha. 1.0 is "
              "LASSO, 0 is Ridge.", false, 0.3, "0.3");
    cmd.add(alpha_arg);

    TCLAP::ValueArg<double>
    flmin_arg("l", "min-lambda", "The minimum lambda as a fraction "
              "of the maximum.", false, 0.3, "0.3");
    cmd.add(flmin_arg);

    TCLAP::ValueArg<unsigned>
    verbosity_arg("v", "verbosity", "Verbosity level (lower is less verbose)",
                  false, 3, "3");
    cmd.add(verbosity_arg);

    TCLAP::SwitchArg switch_scale("s", "scale", "Transform data to z-scores",
                                  cmd, false);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.parse(argc, argv);

    infile = infile_arg.getValue();
    gene_file = genefile_arg.getValue();
    targets_file = targets_arg.getValue();
    bs = batchsize_arg.getValue();
    outfile = outfile_arg.getValue();
    do_scale = switch_scale.getValue();
    ensemble_size = ensemble_arg.getValue();
    min_sample_size = min_exp_arg.getValue();
    max_sample_size = max_exp_arg.getValue();
    predictor_sample_size_min = min_pred_arg.getValue();
    predictor_sample_size_max = max_pred_arg.getValue();
    nlam = nlam_arg.getValue();
    flmin = flmin_arg.getValue();
    alpha = alpha_arg.getValue();
    force = switch_force.getValue();
    verbosity = verbosity_arg.getValue();

    if (targets_file != "") mode = EL_PARTIAL;
    log.set_log_level(verbosity);

  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    log << e.error() << " for arg " << e.argId() << '\n';
    log.log(LOG_ERR);
    return EINVAL;
  }

  // Check all kinds of FS problems that may arise
  if (rank == 0)
  {
    try
    {
      outfile = to_absolute(outfile);
      infile = to_absolute(infile);
      gene_file = to_absolute(gene_file);
      std::string tempdir = dirname(outfile) + "/.seidr_tmp_elnet";

      if (! file_exists(dirname(outfile)) )
        throw std::runtime_error("Directory does not exist: " + dirname(outfile));

      if (dir_exists(tempdir) && force)
      {
        log << "Removing previous temp files.\n";
        log.send(LOG_WARN);
        fs::remove_all(tempdir);
      }

      if (! create_directory(dirname(outfile), ".seidr_tmp_elnet") )
        throw std::runtime_error("Cannot create tmp dir in: " + dirname(outfile));

      if (! file_exists(infile) )
        throw std::runtime_error("File does not exist: " + infile);

      if (! regular_file(infile) )
        throw std::runtime_error("Not a regular file: " + infile);

      if (! file_can_read(infile) )
        throw std::runtime_error("Cannot read: " + infile);

      if (! file_exists(gene_file) )
        throw std::runtime_error("File does not exist: " + gene_file);

      if (! file_can_read(gene_file) )
        throw std::runtime_error("Cannot read: " + gene_file);

      if (mode == EL_PARTIAL)
      {
        targets_file = to_absolute(targets_file);
        if (! file_exists(targets_file) )
          throw std::runtime_error("File does not exist: " + targets_file);

        if (! file_can_read(targets_file))
          throw std::runtime_error("Cannot read: " + targets_file);
      }

      if (! force && file_exists(outfile))
        throw std::runtime_error("File exists: " + outfile);

    }
    catch (std::runtime_error& e)
    {
      log << e.what() << '\n';
      log.log(LOG_ERR);
      return EINVAL;
    }
  }
  // All threads wait until checks are done
  MPI_Barrier(MPI_COMM_WORLD);

  try
  {
    // Get input files
    gene_matrix.load(infile);
    verify_matrix(gene_matrix);
    genes = read_genes(gene_file, row_delim, field_delim);

    if (genes.size() != gene_matrix.n_cols)
      throw std::runtime_error("There must be as many gene names as columns "
                               "in the expression matrix.");

    if (do_scale)
      scale(gene_matrix);

    if (mode == EL_PARTIAL)
      targets = read_genes(targets_file, row_delim, field_delim);

    if (min_sample_size == 0)
      min_sample_size = gene_matrix.n_rows * 0.2;
    if (max_sample_size == 0)
      max_sample_size = gene_matrix.n_rows * 0.8;
    if (predictor_sample_size_min == 0)
      predictor_sample_size_min = (gene_matrix.n_cols - 1) * 0.2;
    if (predictor_sample_size_max == 0)
      predictor_sample_size_max = (gene_matrix.n_cols - 1) * 0.8;

    // Check if sampling settings are sane
    if (min_sample_size > max_sample_size)
      throw std::runtime_error("Minimum experiment sample size can't be "
                               "larger than maximum");
    if (max_sample_size > gene_matrix.n_rows)
      throw std::runtime_error("Maximum experiment sample size can't be "
                               "larger than number of experiments");
    if (predictor_sample_size_min > predictor_sample_size_max)
      throw std::runtime_error("Minimum predictor sample size can't be "
                               "larger than maximum");
    if (predictor_sample_size_max > gene_matrix.n_cols - 1)
      throw std::runtime_error("Maximum predictor sample size can't be "
                               "larger than the number of genes - 1");

    if (min_sample_size == 0 ||
        max_sample_size == 0 ||
        predictor_sample_size_min == 0 ||
        predictor_sample_size_max == 0)
      throw std::runtime_error("None of the sampling settings should be 0");

    switch (mode) {
    case EL_FULL:
      el_full(gene_matrix, genes, bs, outfile,
              min_sample_size, max_sample_size,
              predictor_sample_size_min,
              predictor_sample_size_max,
              ensemble_size, alpha, flmin, nlam);
      break;

    case EL_PARTIAL:
      el_partial(gene_matrix, genes, bs, targets, outfile,
                 min_sample_size, max_sample_size,
                 predictor_sample_size_min,
                 predictor_sample_size_max,
                 ensemble_size, alpha, flmin, nlam);
      break;

    default:
      return 1;
    }
  }
  catch (const std::runtime_error& e)
  {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  }
  catch (const std::exception& e)
  {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  }

  return 0;

}
