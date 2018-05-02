// Seidr
#include <common.h>
#include <genie3-fun.h>
#include <fs.h>
#include <mpims.h>
// External
#include <mpi.h>
#include <cerrno>
#include <string>
#include <vector>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#define GENIE3_FULL 0
#define GENIE3_PARTIAL 1

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
  size_t mode = GENIE3_FULL;
  std::string outfile;
  std::vector<std::string> genes;
  std::vector<std::string> targets;
  seidr_uword_t ntree;
  seidr_uword_t mtry;
  seidr_uword_t min_node_size;
  double alpha;
  double minprop;

  // Define program options
  try
  {
    TCLAP::CmdLine cmd("GENIE3 implementation for Seidr", ' ', version);

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

    TCLAP::ValueArg<unsigned long>
    batchsize_arg("b", "batch-size", "Number of genes in batch", false, 20,
                  "20");
    cmd.add(batchsize_arg);

    TCLAP::ValueArg<std::string>   
    outfile_arg("o", "outfile", "Output file path", false, "edgelist.tsv", 
                "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::ValueArg<unsigned long> 
    ntree_arg("n", "ntree", "Number of random forest trees to grow", false, 1000,
        "1000");
    cmd.add(ntree_arg);

    TCLAP::ValueArg<unsigned long> 
    mtry_arg("m", "mtry", "Number of features to sample in each tree", false, 0,
                                            "sqrt(ngenes)");
    cmd.add(mtry_arg);

    TCLAP::ValueArg<unsigned long> 
    min_node_size_arg("N", "min-node-size", "Minimum node size", false, 0,
        "5");
    cmd.add(min_node_size_arg);

    TCLAP::ValueArg<double> 
    alpha_arg("a", "alpha", "Alpha value for random forests", false, 0.5,
                                      "0.5");
    cmd.add(alpha_arg);

    TCLAP::ValueArg<double> 
    minprop_arg("p", "min-prop", "Minimal proportion in random forest", false, 0.1,
                                        "0.1");
    cmd.add(minprop_arg);

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
    ntree = ntree_arg.getValue();
    mtry = mtry_arg.getValue();
    alpha = alpha_arg.getValue();
    minprop = minprop_arg.getValue();
    force = switch_force.getValue();
    min_node_size = min_node_size_arg.getValue();

    if (targets_file != "") mode = GENIE3_PARTIAL;

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
      std::string tempdir = dirname(outfile) + "/.seidr_tmp_genie3";

      if (! file_exists(dirname(outfile)) )
        throw std::runtime_error("Directory does not exist: " + dirname(outfile));

      if (dir_exists(tempdir) && force)
      {
        log << "Removing previous temp files.\n";
        log.send(LOG_WARN);
        fs::remove_all(tempdir);
      }

      if (! create_directory(dirname(outfile), ".seidr_tmp_genie3") )
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

      if (mode == GENIE3_PARTIAL)
      {
        targets_file = to_absolute(targets_file);
        if (! file_exists(targets_file) )
          throw std::runtime_error("File does not exist: " + targets_file);

        if (! file_can_read(targets_file))
          throw std::runtime_error("Cannot read: " + targets_file);
      }
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
    genes = read_genes(gene_file, row_delim, field_delim);

    if (genes.size() != gene_matrix.n_cols)
      throw std::runtime_error("There must be as many gene names as columns "
                               "in the expression matrix.");

    if (do_scale)
      scale(gene_matrix);

    if (mode == GENIE3_PARTIAL)
      targets = read_genes(targets_file, row_delim, field_delim);

    switch (mode) {
    case GENIE3_FULL:
      genie3_full(gene_matrix, genes, bs, outfile, ntree, mtry, alpha,
                  min_node_size, minprop);
      break;

    case GENIE3_PARTIAL:
      genie3_partial(gene_matrix, genes, bs, targets, outfile, ntree, mtry,
                     min_node_size, alpha, minprop);
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
