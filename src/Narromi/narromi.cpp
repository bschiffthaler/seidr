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
#include <narromi_fun.h>
#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif
// External
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cerrno>
#include <omp.h>
#include <string>
#include <vector>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

int
main(int argc, char** argv)
{

  SEIDR_MPI_INIT();
  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());

  arma::mat gene_matrix;
  std::vector<std::string> genes;
  std::vector<std::string> targets;

  seidr_narromi_param_t param;

  po::options_description umbrella("Narromi implementation for Seidr");
  po::variables_map vm;

  try {
    po::options_description opt("Common Options");
    opt.add_options()("help,h", "Show this help message")(
      "force,f",
      po::bool_switch(&param.force)->default_value(false),
      "Force overwrite if output already exists")(
      "targets,t",
      po::value<std::string>(&param.targets_file),
      "File containing gene names"
      " of genes of interest. The network will only be"
      " calculated using these as the sources of potential connections.")(
      "outfile,o",
      po::value<std::string>(&param.outfile)
        ->default_value("narromi_scores.tsv"),
      "Output file path")("save-resume",
                          po::value<std::string>(&param.cmd_file),
                          "Path to a file that stores job resume info.")(
      "resume-from",
      po::value<std::string>(&param.cmd_file),
      "Try to resume job from this file.")(
      "verbosity,v",
      po::value<unsigned>(&param.verbosity)->default_value(3),
      "Verbosity level (lower is less verbose)");

    po::options_description algopt("NARROMI Options");
    algopt.add_options()(
      "algorithm,m",
      po::value<std::string>(&param.al)->default_value("simplex"),
      "Method for linear programming "
      "optimisaton. One of 'interior-point' or 'simplex'.")(
      "alpha,a",
      po::value<double>(&param.alpha)->default_value(0.05, "0.05"),
      "Alpha cutoff for MI selection.");

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

    po::options_description req("Required");
    req.add_options()("infile,i",
                      po::value<std::string>(&param.infile)->required(),
                      "The expression table (without headers)")(
      "genes,g",
      po::value<std::string>(&param.gene_file)->required(),
      "File containing gene names");

    umbrella.add(req).add(algopt).add(mpiopt).add(ompopt).add(opt);

    po::store(po::command_line_parser(argc, argv).options(umbrella).run(), vm);
  } catch (std::exception& e) {
    log << "Argument exception: " << e.what() << '\n';
    log.send(LOG_ERR);
    return 22;
  }

  if (vm.count("help") || argc == 1) {
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

  param.outfile = to_absolute(param.outfile);
  param.infile = to_absolute(param.infile);
  param.gene_file = to_absolute(param.gene_file);
  if (vm.count("targets")) {
    param.targets_file = to_absolute(param.targets_file);
    param.mode = NARROMI_PARTIAL;
  }

  cp_resume<seidr_narromi_param_t> cp_res(param, CPR_M);
  if (vm.count("resume-from") > 0) {
    assert_exists(param.cmd_file);
    cp_res.load(param, param.cmd_file);
  }

  // Check all kinds of FS problems that may arise
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

        if (vm.count("targets")) {
          assert_exists(param.targets_file);
          assert_can_read(param.targets_file);
        }

        if (!param.force)
          assert_no_overwrite(param.outfile);

        if (!vm.count("tempdir"))
          param.tempdir = tempfile(dirname(param.outfile));
        else
          param.tempdir = tempfile(to_absolute(param.tempdir));
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
        assert_arg_constraint<std::string>({ "simplex", "interior-point" },
                                           param.al);
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
  SEIDR_MPI_BARRIER();

  try {
    assert_in_range<int>(param.nthreads, 1, omp_get_max_threads(), "--threads");

#ifndef NARROMI_USE_CLP
    omp_set_num_threads(1);
    if (param.nthreads > 1) {
      log << "Using GLPK backend\n";
      log.log(LOG_INFO);
      log << "GLPK is not thread safe. Either compile seidr with CLP "
          << "as a backend for Narromi or set you desired number of worker "
          << "threads to be MPI tasks (at the expense of memory usage).\n";
      log.log(LOG_WARN);
    }
#else
    omp_set_num_threads(param.nthreads);
    log << "Using CLP backend\n";
    log.log(LOG_INFO);
#endif
    gene_matrix.load(param.infile);
    genes = read_genes(param.gene_file, param.row_delim, param.field_delim);
    verify_matrix(gene_matrix, genes);
    if (param.mode == NARROMI_PARTIAL) {
      targets =
        read_genes(param.targets_file, param.row_delim, param.field_delim);
    }
    if (vm["batch-size"].defaulted()) {
      if (param.mode == NARROMI_PARTIAL) {
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
  } catch (std::runtime_error& e) {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return errno;
  }

  if (rank == 0) {
    if (vm.count("save-resume") > 0) {
      cp_res.save(param.cmd_file, param);
    }
  }

  try {
    switch (param.mode) {
      case NARROMI_FULL:
        full_narromi(gene_matrix, genes, param);
        break;

      case NARROMI_PARTIAL:
        partial_narromi(gene_matrix, genes, targets, param);
        break;

      default:
        return 1;
    }
  } catch (std::runtime_error& e) {
    log << e.what() << '\n';
    log.log(LOG_ERR);
    return errno;
  }

  return 0;
}
