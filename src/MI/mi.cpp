// Seidr
#include <common.h>
#include <fs.h>
#include <mpims.h>
#include <mi_fun.h>
// External
#include <mpi.h>
#include <cerrno>
#include <string>
#include <vector>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  seidr_mpi_logger log;

  arma::mat gene_matrix;
  std::string infile;
  size_t bs;
  std::string outfile;
  std::string mi_file;
  size_t spline_order;
  size_t num_bins;
  std::string mode;
  bool force = false;
  std::string targets_file;
  std::string genes_file;

  // Define program options
  try
  {
    TCLAP::CmdLine cmd("MI based algorithms for seidr.\n"
                       "These include post processing schemes like CLR "
                       "or DPI (==ARACNE).",
                       ' ', version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<unsigned long>
    batchsize_arg("B", "batch-size", "Number of genes in MPI batch", false, 20,
                  "20");
    cmd.add(batchsize_arg);

    TCLAP::ValueArg<unsigned long>
    spline_arg("s", "spline", "Spline order", false, 3,
               "3");
    cmd.add(spline_arg);

    TCLAP::ValueArg<std::string>
    targets_arg("t", "targets", "File containing gene names"
                " of genes of interest. The network will only be"
                " calculated using these as the sources of potential connections.",
                false, "", "");
    cmd.add(targets_arg);

    TCLAP::ValueArg<std::string>
    genefile_arg("g", "genes", "File containing gene names", true, "",
                 "");
    cmd.add(genefile_arg);

    TCLAP::ValueArg<unsigned long>
    bins_arg("b", "bins", "Number of bins (0 = auto detection)", false, 0,
             "auto");
    cmd.add(bins_arg);

    TCLAP::ValueArg<std::string>
    outfile_arg("o", "outfile", "Output file path", false, "edgelist.tsv",
                "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::ValueArg<std::string>
    mi_file_arg("M", "mi-out", "Save raw MI to this file", false, "", "");
    cmd.add(mi_file_arg);

    std::vector<std::string> modes{"CLR", "ARACNE", "RAW"};
    TCLAP::ValuesConstraint<std::string> modes_constraint(modes);
    TCLAP::ValueArg<std::string>
    mode_arg("m", "mode", "Post processing <RAW>", false, "RAW",
             &modes_constraint);
    cmd.add(mode_arg);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.parse(argc, argv);

    infile = infile_arg.getValue();
    bs = batchsize_arg.getValue();
    outfile = outfile_arg.getValue();
    mi_file = mi_file_arg.getValue();
    spline_order = spline_arg.getValue();
    num_bins = bins_arg.getValue();
    mode = mode_arg.getValue();
    force = switch_force.getValue();
    genes_file = genefile_arg.getValue();
    targets_file = targets_arg.getValue();
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
      if (mi_file != "")
        mi_file = to_absolute(mi_file);
      infile = to_absolute(infile);
      std::string tempdir = dirname(outfile) + "/.seidr_tmp_mi_" + mode;

      if (! file_exists(dirname(outfile)) )
        throw std::runtime_error("Directory does not exist: " + dirname(outfile));

      if (dir_exists(tempdir) && force)
      {
        log << "Removing previous temp files.\n";
        log.send(LOG_WARN);
        fs::remove_all(tempdir);
      }

      if (! create_directory(dirname(outfile), ".seidr_tmp_mi_" + mode) )
        throw std::runtime_error("Cannot create tmp dir in: " + dirname(outfile));

      if (! file_exists(infile) )
        throw std::runtime_error("File does not exist: " + infile);

      if (! regular_file(infile) )
        throw std::runtime_error("Not a regular file: " + infile);

      if (! file_can_read(infile) )
        throw std::runtime_error("Cannot read: " + infile);

      if (! file_can_read(genes_file) )
        throw std::runtime_error("Cannot read: " + genes_file);

      if (! file_can_read(targets_file) && targets_file != "" )
        throw std::runtime_error("Cannot read: " + targets_file);

      if ((! force) && file_exists(outfile))
        throw std::runtime_error("File exists: " + outfile);

      if ((! force) && mi_file != "")
        if (file_exists(mi_file))
          throw std::runtime_error("File exists: " + mi_file);


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
    if (num_bins == 0)
    {
      arma::uvec bc = bin_count(gene_matrix, 1);

      double s = arma::stddev(bc);
      num_bins = arma::median(bc);
      if (num_bins < 2)
        throw std::runtime_error("Autodetected bin count " + std::to_string(num_bins) +
                                 " is too low.");
      else if (num_bins > 15)
      {
        log << "Bin count " << num_bins
            << " is high and might use a high amount of memory.\n";
        log.log(LOG_WARN);
      }
      else
      {
        log << "Detected bin count " << num_bins << ", stddev: " <<
            s << '\n';
        log.log(LOG_INFO);
      }
    }

    if (gene_matrix.n_cols < gene_matrix.n_rows)
    {
      log << "Your expression matrix indicates that you " <<
          "have fewer genes than samples. This is rarely the case. " <<
          "Please make sure genes are columns, samples are rows.\n";
      log.log(LOG_WARN);
    }
    char m = 0;
    if (mode == "CLR")
      m = 0;
    else if (mode == "ARACNE")
      m = 1;
    else if (mode == "RAW")
      m = 2;

    std::vector<std::string> genes;
    genes = read_genes(genes_file);
    std::vector<std::string> targets;
    if (targets_file != "")
      targets = read_genes(targets_file);
    mi_full(gene_matrix, spline_order, num_bins, bs, outfile, m, mi_file,
            genes, targets);
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
