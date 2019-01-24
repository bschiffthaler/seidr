// Seidr
#include <common.h>
#include <fs.h>
#include <mpims.h>
#include <plsnet-fun.h>
// External
#include <mpi.h>
#include <cerrno>
#include <string>
#include <vector>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>

#define PLSNET_FULL 0
#define PLSNET_PARTIAL 1

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
  size_t mode = PLSNET_FULL;
  std::string outfile;
  std::vector<std::string> genes;
  std::vector<std::string> targets;
  arma::uword predictor_sample_size;
  arma::uword ensemble_size;
  arma::uword ncomp;

  // Define program options
  try
  {
    TCLAP::CmdLine cmd("PLSNET implementation for Seidr", ' ',
                       version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<std::string>
    genefile_arg("g", "genes", "File containing gene names (columns in infile)",
                 true, "", "");
    cmd.add(genefile_arg);

    TCLAP::ValueArg<std::string>
    targets_arg("t", "targets", "File containing gene names of genes of interest."
                " The network will only be calculated using these as"
                " the sources of potential connections.",
                false, "", "");
    cmd.add(targets_arg);

    TCLAP::ValueArg<unsigned long>
    batchsize_arg("b", "batch-size", "Number of genes in MPI batch", false, 20,
                  "20");
    cmd.add(batchsize_arg);

    TCLAP::ValueArg<std::string>
    outfile_arg("o", "outfile", "Output file path", false, "edgelist.tsv",
                "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::ValueArg<seidr_uword_t>
    ensemble_arg("e", "ensemble", "The ensemble size", false, 1000, "1000");
    cmd.add(ensemble_arg);

    TCLAP::ValueArg<seidr_uword_t>
    pred_arg("p", "predictor-size", "The number of predictors to be sampled.",
             false, 0, "sqrt(genes)");
    cmd.add(pred_arg);

    TCLAP::ValueArg<seidr_uword_t>
    ncomp_arg("c", "components", "The number of PLS components to be considered.",
              false, 5, "5");
    cmd.add(ncomp_arg);

    TCLAP::SwitchArg
    switch_scale("s", "scale", "Transform data to z-scores", cmd, false);

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
    predictor_sample_size = pred_arg.getValue();
    ncomp = ncomp_arg.getValue();
    force = switch_force.getValue();

    if (targets_file != "") mode = PLSNET_PARTIAL;

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
      std::string tempdir = dirname(outfile) + "/.seidr_tmp_plsnet";

      if (! file_exists(dirname(outfile)) )
        throw std::runtime_error("Directory does not exist: " + dirname(outfile));

      if (dir_exists(tempdir) && force)
      {
        log << "Removing previous temp files.\n";
        log.send(LOG_WARN);
        fs::remove_all(tempdir);
      }

      if (! create_directory(dirname(outfile), ".seidr_tmp_plsnet") )
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

      if (mode == PLSNET_PARTIAL)
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

    if (mode == PLSNET_PARTIAL)
      targets = read_genes(targets_file, row_delim, field_delim);

    if (predictor_sample_size == 0)
    {
      double N = gene_matrix.n_cols - 1;
      N = sqrt(N);
      try
      {
        predictor_sample_size = boost::numeric_cast<arma::uword>(N);
      }
      catch (boost::bad_numeric_cast& e)
      {
        throw std::runtime_error(e.what());
      }
    }

    log << "Sample size: " << predictor_sample_size << '\n';
    log.send(LOG_INFO);

    switch (mode) {
    case PLSNET_FULL:
      plsnet_full(gene_matrix, genes, bs, outfile,
                  predictor_sample_size,
                  ensemble_size, ncomp);
      break;

    case PLSNET_PARTIAL:
      plsnet_partial(gene_matrix, genes, bs, targets, outfile,
                     predictor_sample_size,
                     ensemble_size, ncomp);
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
