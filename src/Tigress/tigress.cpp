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
#include <tiglm.h>
#include <fs.h>
#include <mpiomp.h>
// External
#include <mpi.h>
#include <cerrno>
#include <string>
#include <vector>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  seidr_mpi_logger log;

  arma::mat gene_matrix;
  std::vector<std::string> genes;
  std::vector<std::string> targets;

  seidr_tigress_param_t param;

  po::options_description umbrella("TIGRESS implementation for Seidr");

  po::options_description opt("Common Options");
  opt.add_options()
  ("help,h", "Show this help message")
  ("targets,t", po::value<std::string>(&param.targets_file),
   "File containing gene names"
   " of genes of interest. The network will only be"
   " calculated using these as the sources of potential connections.")
  ("outfile,o",
   po::value<std::string>(&param.outfile)->default_value("elnet_scores.tsv"),
   "Output file path")
  ("verbosity,v",
   po::value<unsigned>(&param.verbosity)->default_value(3),
   "Verbosity level (lower is less verbose)")
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite if output already exists");

  po::options_description algopt("TIGRESS Options");
  algopt.add_options()
  ("scale,s", po::bool_switch(&param.do_scale)->default_value(false),
   "Transform data to z-scores")
  ("nlambda,n",
   po::value<seidr_uword_t>(&param.nsteps)->default_value(10),
   "The maximum number of lambda values")
  ("min-lambda,l",
   po::value<double>(&param.fmin)->default_value(0.3, "0.3"),
   "The minimum lambda as a fraction of the maximum.");
  
  po::options_description mpiopt("MPI Options");
  mpiopt.add_options()
  ("batch-size,B", po::value<uint64_t>(&param.bs)->default_value(20),
   "Number of genes in MPI batch")
  ("tempdir,T",
   po::value<std::string>(&param.tempdir),
   "Temporary directory path");

  po::options_description ompopt("OpenMP Options");
    ompopt.add_options()
    ("threads,O", po::value<int>(&param.nthreads)->
     default_value(omp_get_max_threads()),
     "Number of OpenMP threads per MPI task");

  po::options_description bootopt("Bootstrap Options");
  bootopt.add_options()
  ("ensemble,e",
   po::value<seidr_uword_t>(&param.boot)->default_value(1000),
   "The ensemble size");

  po::options_description req("Required Options");
  req.add_options()
  ("infile,i", po::value<std::string>(&param.infile)->required(),
   "The expression table (without headers)")
  ("genes,g", po::value<std::string>(&param.gene_file)->required(),
   "File containing gene names");

  umbrella.add(req).add(algopt).add(bootopt).add(mpiopt).add(ompopt).add(opt);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(umbrella).run(), vm);

  if (vm.count("help") || argc == 1)
  {
    std::cerr << umbrella << '\n';
    return 1;
  }

  log.set_log_level(5);

  try
  {
    po::notify(vm);
  }
  catch (std::exception& e)
  {
    log << "Argument exception: " << e.what() << '\n';
    log.send(LOG_ERR);
    return 22;
  }

  log.set_log_level(param.verbosity);

  if (vm.count("targets"))
    param.mode = TIGLM_PARTIAL;

  // Normalize paths
  param.outfile = to_absolute(param.outfile);
  param.infile = to_absolute(param.infile);
  param.gene_file = to_absolute(param.gene_file);

  if (param.mode == TIGLM_PARTIAL)
  {
    param.targets_file = to_absolute(param.targets_file);
  }

  // Check all kinds of FS problems that may arise in the master thread
  if (rank == 0)
  {
    try
    {
      assert_exists(dirname(param.outfile));
      assert_exists(param.infile);
      assert_is_regular_file(param.infile);
      assert_exists(param.gene_file);
      assert_can_read(param.gene_file);
      assert_can_read(param.infile);

      if (param.mode == TIGLM_PARTIAL)
      {
        assert_exists(param.targets_file);
        assert_can_read(param.targets_file);
      }

      if (! param.force)
        assert_no_overwrite(param.outfile);

      if (! vm.count("tempdir"))
        param.tempdir = tempfile(dirname(param.outfile));
      else
        param.tempdir = tempfile(to_absolute(param.tempdir));
      if (dir_exists(param.tempdir))
      {
        if (param.force)
        {
          log << "Removing previous temp files.\n";
          log.log(LOG_WARN);
          fs::remove_all(param.tempdir);
        }
        else
        {
          throw std::runtime_error("Dir exists: " + param.tempdir);
        }
      }
      else
      {
        create_directory(param.tempdir);
      }

      assert_dir_is_writeable(param.tempdir);
      mpi_sync_tempdir(&param.tempdir);
    }
    catch (std::runtime_error& e)
    {
      log << e.what() << '\n';
      log.log(LOG_ERR);
      return EINVAL;
    }
  }
  else
  {
    mpi_sync_tempdir(&param.tempdir);
  }
  
  // All threads wait until checks are done
  MPI_Barrier(MPI_COMM_WORLD);

  try
  {
    assert_in_range<int>(param.nthreads, 1, omp_get_max_threads(),
                         "--threads");
    omp_set_num_threads(param.nthreads);
    // Get input files
    gene_matrix.load(param.infile);
    genes = read_genes(param.gene_file, param.row_delim, param.field_delim);

    verify_matrix(gene_matrix, genes);

    if (param.do_scale)
      scale(gene_matrix);

    if (param.mode == TIGLM_PARTIAL)
      targets = read_genes(param.targets_file, param.row_delim, param.field_delim);

    switch (param.mode) {
    case TIGLM_FULL:
      tiglm_full(gene_matrix, genes, param);
      break;

    case TIGLM_PARTIAL:
      tiglm_partial(gene_matrix, genes, targets, param);
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
