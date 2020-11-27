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
#include <tiglm.h>
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
#include <string>
#include <vector>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

int
main(int argc, char** argv)
{

  SEIDR_MPI_INIT();

  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());

  try {

    arma::mat gene_matrix;
    std::vector<std::string> genes;
    std::vector<std::string> targets;

    seidr_tigress_param_t param;

    po::options_description umbrella("TIGRESS implementation for Seidr");
    po::variables_map vm;

    po::options_description opt("Common Options");
    opt.add_options()("help,h", "Show this help message")(
      "targets,t",
      po::value<std::string>(&param.targets_file),
      "File containing gene names"
      " of genes of interest. The network will only be"
      " calculated using these as the sources of potential connections.")(
      "outfile,o",
      po::value<std::string>(&param.outfile)
        ->default_value("tigress_scores.tsv"),
      "Output file path")("save-resume",
                          po::value<std::string>(&param.cmd_file),
                          "Path to a file that stores job resume info.")(
      "resume-from",
      po::value<std::string>(&param.cmd_file),
      "Try to resume job from this file.")(
      "verbosity,v",
      po::value<unsigned>(&param.verbosity)
        ->default_value(TIGRESS_DEF_VERBOSITY),
      "Verbosity level (lower is less verbose)")(
      "force,f",
      po::bool_switch(&param.force)->default_value(false),
      "Force overwrite if output already exists")(
      "version,V", "Print the program version");

    po::options_description algopt("TIGRESS Options");
    algopt.add_options()(
      "scale,s", "(deprecated) Transform data to z-scores")(
      "no-scale,S",
      po::bool_switch(&param.do_scale)->default_value(true),
      "Do not transform data to z-scores")(
      "nlambda,n",
      po::value<seidr_uword_t>(&param.nsteps)
        ->default_value(TIGRESS_DEF_NSTEPS),
      "The maximum number of lambda values")(
      "min-lambda,l",
      po::value<double>(&param.fmin)
        ->default_value(TIGRESS_DEF_FMIN, _XSTR(TIGRESS_DEF_FMIN)),
      "The minimum lambda as a fraction of the maximum.")(
      "allow-low-var",
      po::bool_switch(&param.allow_low_var)->default_value(false),
      "Disable low variance filter when subsampling data (at your own risk)");

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
    bootopt.add_options()("ensemble,e",
                          po::value<seidr_uword_t>(&param.boot)
                            ->default_value(TIGRESS_DEF_ENSEMBLE),
                          "The ensemble size");

    po::options_description req("Required Options");
    req.add_options()("infile,i",
                      po::value<std::string>(&param.infile)->required(),
                      "The expression table (without headers)")(
      "genes,g",
      po::value<std::string>(&param.gene_file)->required(),
      "File containing gene names");

    umbrella.add(req).add(algopt).add(bootopt).add(mpiopt).add(ompopt).add(opt);

    po::store(po::command_line_parser(argc, argv).options(umbrella).run(), vm);

    if (vm.count("help") != 0 || argc == 1) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    if (vm.count("version") > 0) {
      std::cout << _XSTR(SEIDR_VERSION) << '\n';
      SEIDR_MPI_FINALIZE();
      return 0;
    }

    seidr_mpi_logger::set_log_level(LOG_DEBUG);

    po::notify(vm);

    seidr_mpi_logger::set_log_level(param.verbosity);

    if (vm.count("targets") != 0) {
      param.mode = TIGLM_PARTIAL;
    }

    if (vm.count("scale") > 0) {
      log << "--scale is deprecated as it is now default. Use --no-scale"
          << " to turn scaling off.\n";
      log.log(LOG_WARN);
    }

    // Normalize paths
    param.outfile = to_absolute(param.outfile);
    param.infile = to_absolute(param.infile);
    param.gene_file = to_absolute(param.gene_file);

    if (param.mode == TIGLM_PARTIAL) {
      param.targets_file = to_absolute(param.targets_file);
    }

    cp_resume<seidr_tigress_param_t> cp_res(param, CPR_M);
    if (vm.count("resume-from") > 0) {
      assert_exists(param.cmd_file);
      cp_res.load(param, param.cmd_file);
    }

    // Check all kinds of FS problems that may arise in the master thread
    if (rank == 0) {
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
        assert_no_cr(param.gene_file);
        assert_no_cr(param.infile);

        if (param.mode == TIGLM_PARTIAL) {
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
    } else {
      mpi_sync_tempdir(&param.tempdir);
      if (vm.count("resume-from") > 0) {
        mpi_sync_cpr_vector(&param.good_idx);
      }
    }

    // All threads wait until checks are done
    SEIDR_MPI_BARRIER();

    assert_in_range<int>(param.nthreads, 1, omp_get_max_threads(), "--threads");
    log << "Setting " << param.nthreads << " OMP threads\n";
    log.send(LOG_DEBUG);
    omp_set_num_threads(param.nthreads);

    // Get input files
    log << "Loading " << param.infile << '\n';
    log.send(LOG_DEBUG);
    gene_matrix.load(param.infile);

    log << "Reading " << param.gene_file << '\n';
    log.send(LOG_DEBUG);
    genes = read_genes(param.gene_file, param.row_delim, param.field_delim);

    log << "Checking if matrix looks OK\n";
    log.send(LOG_DEBUG);
    verify_matrix(gene_matrix, genes);

    if (param.do_scale) {
      scale(gene_matrix);
    }

    if (param.mode == TIGLM_PARTIAL) {
      targets =
        read_genes(param.targets_file, param.row_delim, param.field_delim);
    }

    if (vm["batch-size"].defaulted()) {
      if (param.mode == TIGLM_PARTIAL) {
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

    if (rank == 0) {
      if (vm.count("save-resume") > 0) {
        cp_res.save(param.cmd_file, param);
      }
    }

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
  } catch (const po::error& e) {
    log << "[Argument Error]: " << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  } catch (const std::runtime_error& e) {
    log << "[Runtime Error]: " << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  } catch (const std::exception& e) {
    log << "[Generic Error]: " << e.what() << '\n';
    log.log(LOG_ERR);
    return 1;
  }

  return 0;
}
