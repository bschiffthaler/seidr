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
#include <cp_resume.h>
#include <fs.h>
#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif
#include <elnet-fun.h>
// External
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <cerrno>
#include <omp.h>
#include <string>
#include <vector>

namespace fs = boost::filesystem;
namespace po = boost::program_options;
using boost::numeric_cast;

int
main(int argc, char** argv)
{

  SEIDR_MPI_INIT();

  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());

  seidr_elnet_param_t param;
  po::variables_map vm;
  po::options_description umbrella("Elastic net ensemble implementation for"
                                   " Seidr");
  try {
    po::options_description opt("Common Options");
    opt.add_options()("help,h", "Show this help message")(
      "targets,t",
      po::value<std::string>(&param.targets_file),
      "File containing gene names"
      " of genes of interest. The network will only be"
      " calculated using these as the sources of potential connections.")(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("elnet_scores.tsv"),
      "Output file path")("save-resume",
                          po::value<std::string>(&param.cmd_file),
                          "Path to a file that stores job resume info.")(
      "resume-from",
      po::value<std::string>(&param.cmd_file),
      "Try to resume job from this file.")(
      "verbosity,v",
      po::value<unsigned>(&param.verbosity)->default_value(3),
      "Verbosity level (lower is less verbose)")(
      "force,f",
      po::bool_switch(&param.force)->default_value(false),
      "Force overwrite if output already exists");

    po::options_description algopt("ElNet Options");
    algopt.add_options()("scale,s",
                         po::bool_switch(&param.do_scale)->default_value(false),
                         "Transform data to z-scores")(
      "nlambda,n",
      po::value<seidr_uword_t>(&param.nlam)->default_value(10),
      "The maximum number of lambda values")(
      "min-lambda,l",
      po::value<double>(&param.flmin)->default_value(0.3, "0.3"),
      "The minimum lambda as a fraction of the maximum.")(
      "alpha,a",
      po::value<double>(&param.alpha)->default_value(0.3, "0.3"),
      "The elastic net mixing value alpha. 1.0 is "
      "LASSO, 0 is Ridge.");

    po::options_description mpiopt("MPI Options");
    mpiopt.add_options()("batch-size,B",
                         po::value<uint64_t>(&param.bs)->default_value(0),
                         "Number of genes in MPI batch")(
      "tempdir,T",
      po::value<std::string>(&param.tempdir),
      "Temporary directory path");

    po::options_description ompopt("OpenMP Options");
    ompopt.add_options()(
      "threads,O",
      po::value<int>(&param.nthreads)->default_value(omp_get_max_threads()),
      "Number of OpenMP threads per MPI task");

    po::options_description bootopt("Bootstrap Options");
    bootopt.add_options()(
      "ensemble,e",
      po::value<seidr_uword_t>(&param.ensemble_size)->default_value(1000),
      "The ensemble size")(
      "min-predictor-size,p",
      po::value<seidr_uword_t>(&param.predictor_sample_size_min)
        ->default_value(0, "20% of predictors"),
      "The minimum number of predictors to be sampled.")(
      "max-predictor-size,P",
      po::value<seidr_uword_t>(&param.predictor_sample_size_max)
        ->default_value(0, "80% of predictors"),
      "The maximum number of predictors to be sampled")(
      "min-experiment-size,x",
      po::value<seidr_uword_t>(&param.min_sample_size)
        ->default_value(0, "20% of experiments"),
      "The minimum number of experiments to be sampled")(
      "max-experiment-size,X",
      po::value<seidr_uword_t>(&param.max_sample_size)
        ->default_value(0, "80% of experiments"),
      "The maximum number of experiments to be sampled");

    po::options_description req("Required Options");
    req.add_options()("infile,i",
                      po::value<std::string>(&param.infile)->required(),
                      "The expression table (without headers)")(
      "genes,g",
      po::value<std::string>(&param.gene_file)->required(),
      "File containing gene names");

    umbrella.add(req).add(algopt).add(bootopt).add(mpiopt).add(ompopt).add(opt);

    po::store(po::command_line_parser(argc, argv).options(umbrella).run(), vm);
  } catch (std::exception& e) {
    log << "Argument exception: " << e.what() << '\n';
    log.send(LOG_ERR);
    return 22;
  }

  if (vm.count("help") != 0 || argc == 1) {
    std::cerr << umbrella << '\n';
    return 1;
  }

  log.set_log_level(5);

  try {
    po::notify(vm);
  } catch (std::exception& e) {
    log << "Argument exception: " << e.what() << '\n';
    log.send(LOG_ERR);
    return 22;
  }

  log.set_log_level(param.verbosity);

  if (vm.count("targets") != 0) {
    param.mode = EL_PARTIAL;
  }

  // Normalize paths
  param.outfile = to_absolute(param.outfile);
  param.infile = to_absolute(param.infile);
  param.gene_file = to_absolute(param.gene_file);

  if (param.mode == EL_PARTIAL) {
    param.targets_file = to_absolute(param.targets_file);
  }

  cp_resume<seidr_elnet_param_t> cp_res(param, CPR_M);
  if (vm.count("resume-from") > 0) {
    assert_exists(param.cmd_file);
    cp_res.load(param, param.cmd_file);
  }

  // Check all kinds of FS problems that may arise in the master thread
  if (rank == 0) {
    try {
      if (vm.count("resume-from") > 0 && file_exists(param.cmd_file)) {
        log << "Trying to resume from " << param.cmd_file << '\n';
        log.log(LOG_INFO);
        cp_res.resume();
        mpi_sync_tempdir(&param.tempdir);
        // Only need to sync CPR when resuming
        param.good_idx = cp_res.get_good_idx();
        mpi_sync_cpr_vector(&param.good_idx);
      } else {
        assert_exists(dirname(param.outfile));
        assert_exists(param.infile);
        assert_is_regular_file(param.infile);
        assert_exists(param.gene_file);
        assert_can_read(param.gene_file);
        assert_can_read(param.infile);

        if (param.mode == EL_PARTIAL) {
          assert_exists(param.targets_file);
          assert_can_read(param.targets_file);
        }

        if (!param.force) {
          assert_no_overwrite(param.outfile);
        }

        if (vm.count("tempdir") == 0) {
          param.tempdir = tempfile(dirname(param.outfile));
        } else {
          param.tempdir = tempfile(to_absolute(param.tempdir));
        }
        if (dir_exists(param.tempdir)) {
          if (param.force) {
            log << "Removing previous temp files.\n";
            log.log(LOG_WARN);
            fs::remove_all(param.tempdir);
          } else {
            throw std::runtime_error("Dir exists: " + param.tempdir);
          }
        } else {
          create_directory(param.tempdir);
        }

        assert_dir_is_writeable(param.tempdir);
        mpi_sync_tempdir(&param.tempdir);
      }
    } catch (std::runtime_error& e) {
      log << e.what() << '\n';
      log.log(LOG_ERR);
      return EINVAL;
    }
  } else {
    mpi_sync_tempdir(&param.tempdir);
    if (vm.count("resume-from") > 0) {
      mpi_sync_cpr_vector(&param.good_idx);
    }
  }
  // All threads wait until checks are done
  SEIDR_MPI_BARRIER(); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  try {
    assert_in_range<int>(param.nthreads, 1, omp_get_max_threads(), "--threads");
    omp_set_num_threads(param.nthreads);
    arma::mat gene_matrix;
    std::vector<std::string> genes;
    std::vector<std::string> targets;

    // Get input files
    gene_matrix.load(param.infile);
    genes = read_genes(param.gene_file, param.row_delim, param.field_delim);
    verify_matrix(gene_matrix, genes);

    if (param.do_scale) {
      scale(gene_matrix);
    }

    if (param.mode == EL_PARTIAL) {
      targets =
        read_genes(param.targets_file, param.row_delim, param.field_delim);
    }

    if (vm["batch-size"].defaulted()) {
      if (param.mode == EL_PARTIAL) {
        uint64_t s = param.resuming ? (targets.size() - param.good_idx.size())
                                    : targets.size();
        param.bs = guess_batch_size(s, get_mpi_nthread());
      } else {
        uint64_t s = param.resuming ? (genes.size() - param.good_idx.size())
                                    : genes.size();
        param.bs = guess_batch_size(s, get_mpi_nthread());
      }
      log << "Setting batch size to " << param.bs << '\n';
      log.log(LOG_INFO);
    }

    if (param.min_sample_size == 0) {
      param.min_sample_size =
        numeric_cast<arma::uword>(gene_matrix.n_rows * 0.2);
    }
    if (param.max_sample_size == 0) {
      param.max_sample_size =
        numeric_cast<arma::uword>(gene_matrix.n_rows * 0.8);
    }
    if (param.predictor_sample_size_min == 0) {
      param.predictor_sample_size_min =
        numeric_cast<arma::uword>((gene_matrix.n_cols - 1) * 0.2);
    }
    if (param.predictor_sample_size_max == 0) {
      param.predictor_sample_size_max =
        numeric_cast<arma::uword>((gene_matrix.n_cols - 1) * 0.8);
    }

    // Check if sampling settings are sane
    if (param.min_sample_size > param.max_sample_size) {
      throw std::runtime_error("Minimum experiment sample size can't be "
                               "larger than maximum");
    }
    if (param.max_sample_size > gene_matrix.n_rows) {
      throw std::runtime_error("Maximum experiment sample size can't be "
                               "larger than number of experiments");
    }
    if (param.predictor_sample_size_min > param.predictor_sample_size_max) {
      throw std::runtime_error("Minimum predictor sample size can't be "
                               "larger than maximum");
    }
    if (param.predictor_sample_size_max > gene_matrix.n_cols - 1) {
      throw std::runtime_error("Maximum predictor sample size can't be "
                               "larger than the number of genes - 1");
    }

    if (param.min_sample_size == 0 || param.max_sample_size == 0 ||
        param.predictor_sample_size_min == 0 ||
        param.predictor_sample_size_max == 0) {
      throw std::runtime_error("None of the sampling settings should be 0");
    }

    if (rank == 0) {
      if (vm.count("save-resume") > 0) {
        cp_res.save(param.cmd_file, param);
      }
    }

    switch (param.mode) {
      case EL_FULL:
        el_full(gene_matrix, genes, param);
        break;

      case EL_PARTIAL:
        el_partial(gene_matrix, genes, targets, param);
        break;

      default:
        return 1;
    }
  } catch (const std::runtime_error& e) {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  } catch (const std::exception& e) {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  }

  return 0;
}
