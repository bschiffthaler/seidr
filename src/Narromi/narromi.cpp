// Seidr
#include <narromi_fun.h>
#include <common.h>
#include <fs.h>
#include <mpims.h>
// External
#include <mpi.h>
#include <cerrno>
#include <vector>
#include <string>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#define NARROMI_FULL 0
#define NARROMI_PARTIAL 1

namespace fs = boost::filesystem;

int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  seidr_mpi_logger log;

  arma::mat gene_matrix;
  std::string al;
  double alpha = 0.05;
  double t = 0.6;
  bool force = false;
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  size_t mode = NARROMI_FULL;
  char row_delim = '\n';
  char field_delim = '\t';
  size_t batchsize;
  std::vector<std::string> genes;
  std::vector<std::string> targets;
  std::string outfile;

  try
  {
    TCLAP::CmdLine cmd("Narromi implementation for Seidr", ' ', version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<std::string>
    genefile_arg("g", "genes", "File containing gene names (columns in infile)",
                 true, "", "");
    cmd.add(genefile_arg);

    TCLAP::ValueArg<std::string>
    targets_arg("t", "targets", "File containing gene names of genes of interes. "
                "The network will only be calculated using these as the sources of "
                "potential connections.", false, "", "");
    cmd.add(targets_arg);

    TCLAP::ValueArg<std::string>
    outfile_arg("o", "outfile", "Outfile name", false, "edgelist.tsv",
                "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::ValueArg<size_t>
    bs_arg("b", "batch-size", "Batch size for MPI threads", false, 20,
           "20");
    cmd.add(bs_arg);

    TCLAP::ValueArg<std::string>
    algorithm_arg("m", "algorithm", "Method for linear programming "
                  "optimisaton. One of 'interior-point' or 'simplex'.", false,
                  "simplex", "simplex");
    cmd.add(algorithm_arg);

    TCLAP::ValueArg<double> alpha_arg("a", "alpha",
                                      "Alpha cutoff for MI selection (== beta)",
                                      false, 0.05, "0.05");
    cmd.add(alpha_arg);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.parse(argc, argv);

    infile = infile_arg.getValue();
    gene_file = genefile_arg.getValue();
    targets_file = targets_arg.getValue();
    batchsize = bs_arg.getValue();
    al = algorithm_arg.getValue();
    alpha =  alpha_arg.getValue();
    outfile = outfile_arg.getValue();
    force = switch_force.getValue();

    if (targets_file != "") mode = NARROMI_PARTIAL;

  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "[Runtime error]: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  try
  {
    if (rank == 0)
    {
      infile = to_absolute(infile);
      outfile = to_absolute(outfile);
      gene_file = to_absolute(gene_file);
      std::string tempdir = dirname(outfile) + "/.seidr_tmp_narromi";

      if (! file_exists(dirname(outfile)))
        throw std::runtime_error("Directory does not exist: " + dirname(outfile));

      if (dir_exists(tempdir) && force)
      {
        log << "Removing previous temp files.\n";
        log.send(LOG_WARN);
        fs::remove_all(tempdir);
      }

      if (! create_directory(dirname(outfile), ".seidr_tmp_narromi") )
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

      if (mode == NARROMI_PARTIAL)
      {
        targets_file = to_absolute(targets_file);

        if (! file_exists(targets_file) )
          throw std::runtime_error("File does not exist:  " + targets_file);

        if (! file_can_read(targets_file))
          throw std::runtime_error("Cannot read:  " + targets_file);

      }

      if (! force && file_exists(outfile))
        throw std::runtime_error("File exists: " + outfile);

    }
  }
  catch (std::runtime_error& e)
  {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return errno;
  }

  try
  {
    gene_matrix.load(infile);
    genes = read_genes(gene_file, row_delim, field_delim);
    if (mode == NARROMI_PARTIAL)
      targets = read_genes(targets_file, row_delim, field_delim);
    if (genes.size() != gene_matrix.n_cols)
      throw std::runtime_error("There must be as many gene names as columns "
                               "in the expression matrix.");
  }
  catch (std::runtime_error& e)
  {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return errno;
  }

  try
  {
    switch (mode) {
    case NARROMI_FULL:
      full_narromi(gene_matrix, al, alpha, t, genes, batchsize, outfile);
      break;

    case NARROMI_PARTIAL:
      partial_narromi(gene_matrix, al, alpha, t, genes, batchsize, targets, outfile);
      break;

    default:
      return EINVAL;
    }
  }
  catch (std::runtime_error& e)
  {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return errno;
  }

  return 0;
}
