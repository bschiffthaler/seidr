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
#include <mpiomp.h>
#include <svm-fun.h>
// External
#include <mpi.h>
#include <svm.h>
#include <cerrno>
#include <string>
#include <vector>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
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

  seidr_svm_param_t param;

  po::options_description umbrella("NIMEFI SVM implementation for Seidr");

  po::options_description opt("Common Options");
  opt.add_options()
  ("help,h", "Show this help message")
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite if output already exists")
  ("targets,t", po::value<std::string>(&param.targets_file),
   "File containing gene names"
   " of genes of interest. The network will only be"
   " calculated using these as the sources of potential connections.")
  ("outfile,o",
   po::value<std::string>(&param.outfile)->default_value("elnet_scores.tsv"),
   "Output file path")
  ("verbosity,v",
   po::value<unsigned>(&param.verbosity)->default_value(3),
   "Verbosity level (lower is less verbose)");

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
   po::value<seidr_uword_t>(&param.ensemble_size)->default_value(1000),
   "The ensemble size")
  ("min-predictor-size,p",
   po::value<seidr_uword_t>(&param.predictor_sample_size_min)->
   default_value(0, "20% of predictors"),
   "The minimum number of predictors to be sampled.")
  ("max-predictor-size,P",
   po::value<seidr_uword_t>(&param.predictor_sample_size_max)->
   default_value(0, "80% of predictors"),
   "The maximum number of predictors to be sampled")
  ("min-experiment-size,x",
   po::value<seidr_uword_t>(&param.min_sample_size)->
   default_value(0, "20% of experiments"),
   "The minimum number of experiments to be sampled")
  ("max-experiment-size,X",
   po::value<seidr_uword_t>(&param.max_sample_size)->
   default_value(0, "80% of experiments"),
   "The maximum number of experiments to be sampled");

  po::options_description req("Required");
  req.add_options()
  ("infile,i", po::value<std::string>(&param.infile)->required(),
   "The expression table (without headers)")
  ("genes,g", po::value<std::string>(&param.gene_file)->required(),
   "File containing gene names");

  po::options_description svmopt("SVM options");
  svmopt.add_options()
  ("scale,s", po::bool_switch(&param.do_scale)->default_value(false),
   "Transform data to z-scores")
  ("type,y",
   po::value<std::string>(&param.svm_type)->default_value("EPSILON_SVR"),
   "SVM type [NU_SVR, EPSILON_SVR]")
  ("kernel,k",
   po::value<std::string>(&param.kernel)->default_value("LINEAR"),
   "SVM kernel [LINEAR, POLY, RBF, SIGMOID]")
  ("degree,d",
   po::value<int>(&param.svparam.degree)->default_value(3),
   "Polynomial degree (for POLY kernel)")
  ("gamma,G",
   po::value<double>(&param.svparam.gamma)->default_value(0.01, "0.01"),
   "Kernel coefficient for POLY/RBF/SIGMOID kernels")
  ("coef,c",
   po::value<double>(&param.svparam.coef0)->default_value(0.01, "0.01"),
   "Independent term in kernel function (for POLY/SIGMOID kernels)")
  ("nu,n",
   po::value<double>(&param.svparam.nu)->default_value(0.5, "0.5"),
   "nu value (for NU_SVR)")
  ("penalty,C",
   po::value<double>(&param.svparam.C)->default_value(1, "1"),
   "Penalty C value")
  ("tol,l",
   po::value<double>(&param.svparam.eps)->default_value(0.001, "0.001"),
   "Epsilon/tolerance (stopping criterion)")
  ("eps,E",
   po::value<double>(&param.svparam.p)->default_value(0.1, "0.1"),
   "Epsilon (for EPSILON_SVR)")
  ("shrinking,S",
   po::value<int>(&param.svparam.shrinking)->default_value(1),
   "Whether to use the shrinking heuristic [0: off, 1: on]");

  // param.svparam.shrinking = 1;


  umbrella.add(req).add(svmopt).add(bootopt).add(mpiopt).add(ompopt).add(opt);

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
    param.mode = SVM_PARTIAL;

  // Normalize paths
  param.outfile = to_absolute(param.outfile);
  param.infile = to_absolute(param.infile);
  param.gene_file = to_absolute(param.gene_file);

  if (param.mode == SVM_PARTIAL)
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

      if (param.mode == SVM_PARTIAL)
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
      assert_arg_constraint<std::string>({"NU_SVR", "EPSILON_SVR"},
                                         param.svm_type);
      assert_arg_constraint<std::string>({"LINEAR", "POLY", "RBF", "SIGMOID"},
                                         param.kernel);
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

    if (param.mode == SVM_PARTIAL)
      targets = read_genes(param.targets_file, param.row_delim,
                           param.field_delim);

    if (param.min_sample_size == 0)
      param.min_sample_size = gene_matrix.n_rows / 5;
    if (param.max_sample_size == 0)
      param.max_sample_size = 4 * (gene_matrix.n_rows / 5);
    if (param.predictor_sample_size_min == 0)
      param.predictor_sample_size_min = (gene_matrix.n_cols - 1) / 5;
    if (param.predictor_sample_size_max == 0)
      param.predictor_sample_size_max = 4 * ((gene_matrix.n_cols - 1) / 5);

    // Check if sampling settings are sane
    if (param.min_sample_size > param.max_sample_size)
      throw std::runtime_error("Minimum experiment sample size can't be "
                               "larger than maximum");
    if (param.max_sample_size > gene_matrix.n_rows)
      throw std::runtime_error("Maximum experiment sample size can't be "
                               "larger than number of experiments");
    if (param.predictor_sample_size_min > param.predictor_sample_size_max)
      throw std::runtime_error("Minimum predictor sample size can't be "
                               "larger than maximum");
    if (param.predictor_sample_size_max > gene_matrix.n_cols - 1)
      throw std::runtime_error("Maximum predictor sample size can't be "
                               "larger than the number of genes - 1");

    if (param.min_sample_size == 0 ||
        param.max_sample_size == 0 ||
        param.predictor_sample_size_min == 0 ||
        param.predictor_sample_size_max == 0)
      throw std::runtime_error("None of the sampling settings should be 0");

    if (param.kernel == "LINEAR")
    {
      param.svparam.kernel_type = LINEAR;
    }
    else if (param.kernel == "POLY")
    {
      param.svparam.kernel_type = POLY;
    }
    else if (param.kernel == "RBF")
    {
      param.svparam.kernel_type = RBF;
    }
    else if (param.kernel == "SIGMOID")
    {
      param.svparam.kernel_type = SIGMOID;
    }
    else
    {
      throw std::runtime_error("--kernel needs to be one of "
                               "[LINEAR, POLY, RBF, SIGMOID]");
    }

    if (param.svm_type == "EPSILON_SVR")
    {
      param.svparam.svm_type = EPSILON_SVR;
    }
    else if (param.svm_type == "NU_SVR")
    {
      param.svparam.svm_type = NU_SVR;
    }

    param.svparam.cache_size = 256;
    param.svparam.probability = 0;
    param.svparam.nr_weight = 0;
    param.svparam.weight_label = NULL;
    param.svparam.weight = NULL;


    switch (param.mode)
    {
    case SVM_FULL:
      svm_full(gene_matrix, genes, param);
      break;

    case SVM_PARTIAL:
      svm_partial(gene_matrix, genes, targets, param);
      break;

    default:
      return EINVAL;
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
