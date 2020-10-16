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
#include <mi_fun.h>
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

  try {
    arma::mat gene_matrix;

    seidr_mi_param_t param;

    po::options_description umbrella(
      "MI based algorithms for seidr.\n"
      "These include post processing schemes like CLR "
      "or DPI (==ARACNE).");

    po::options_description opt("Common Options");
    po::variables_map vm;

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
      po::value<std::string>(&param.outfile)->default_value("mi_scores.tsv"),
      "Output file path")("save-resume",
                          po::value<std::string>(&param.cmd_file),
                          "Path to a file that stores job resume info.")(
      "resume-from",
      po::value<std::string>(&param.cmd_file),
      "Try to resume job from this file.")(
      "verbosity,v",
      po::value<unsigned>(&param.verbosity)->default_value(MI_DEF_VERBOSITY),
      "Verbosity level (lower is less verbose)")(
      "version,V", po::bool_switch(), "Print the program version");

    po::options_description algopt("MI Options");
    algopt.add_options()("spline,s",
                         po::value<uint64_t>(&param.spline_order)
                           ->default_value(MI_DEF_SPLINE),
                         "Spline order")(
      "bins,b",
      po::value<uint64_t>(&param.num_bins)
        ->default_value(MI_DEF_NUM_BINS, "auto"),
      "Number of bins")("mi-file,M",
                        po::value<std::string>(&param.mi_file),
                        "Save/load raw MI to/from this file")(
      "mode,m",
      po::value<std::string>(&param.mode)->default_value("RAW"),
      "Post processing [RAW, CLR, ARACNE]");

    po::options_description mpiopt("MPI Options");
    mpiopt.add_options()(
      "batch-size,B",
      po::value<uint64_t>(&param.bs)->default_value(MI_DEF_BS),
      "Number of genes in MPI batch")("tempdir,T",
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

    if (vm.count("help") != 0 || argc == 1) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    if (vm.count("version") > 0) {
      std::cout << _XSTR(SEIDR_VERSION) << '\n';
      return 0;
    }

    po::notify(vm);

    seidr_mpi_logger::set_log_level(param.verbosity);

    param.outfile = to_absolute(param.outfile);
    param.infile = to_absolute(param.infile);
    param.gene_file = to_absolute(param.gene_file);
    if (vm.count("mi-file") != 0 && exists(param.mi_file)) {
      param.use_existing = true;
      param.mi_file = to_absolute(param.mi_file);
    }
    if (vm.count("targets") != 0) {
      param.targets_file = to_absolute(param.targets_file);
    }

    cp_resume<seidr_mi_param_t> cp_res(param, CPR_LM);
    if (vm.count("resume-from") != 0) {
      assert_exists(param.cmd_file);
      cp_res.load(param, param.cmd_file);
    }

    // Check all kinds of FS problems that may arise
    if (rank == 0) {

      if (vm.count("resume-from") != 0 && file_exists(param.cmd_file)) {
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

        if (vm.count("targets") != 0) {
          log << "Targets selected, but this algorithm will still need to "
              << "calculate (or load) the full MI matrix first.\n";
          assert_exists(param.targets_file);
          assert_can_read(param.targets_file);
        }

        if (vm.count("mi-file") != 0 && exists(param.mi_file)) {
          assert_can_read(param.mi_file);
        } else if (vm.count("mi-file") != 0) {
          assert_dir_is_writeable(dirname(param.mi_file));
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
        assert_arg_constraint<std::string>({ "RAW", "CLR", "ARACNE" },
                                           param.mode);
        assert_dir_is_writeable(param.tempdir);
        mpi_sync_tempdir(&param.tempdir);
      }

    } else {
      mpi_sync_tempdir(&param.tempdir);
      if (vm.count("resume-from") != 0) {
        mpi_sync_cpr_vector(&param.good_idx);
      }
    }
    // All threads wait until checks are done
    SEIDR_MPI_BARRIER();

    assert_in_range<int>(param.nthreads, 1, omp_get_max_threads(), "--threads");
    omp_set_num_threads(param.nthreads);
    // Get input files
    gene_matrix.load(param.infile);
    std::vector<std::string> genes = read_genes(param.gene_file);
    verify_matrix(gene_matrix, genes);
    if (param.num_bins == 0) {
      auto bc = arma::conv_to<arma::vec>::from(bin_count(gene_matrix, 1));

      double s = arma::stddev(bc);
      param.num_bins = arma::median(bc);
      if (param.num_bins < MI_MIN_BINS) {
        throw std::runtime_error("Autodetected bin count " +
                                 std::to_string(param.num_bins) +
                                 " is too low.");
      }
      if (param.num_bins > MI_WARN_BINS) {
        log << "Bin count " << param.num_bins << ", stddev: " << s
            << " is high and might use a high amount of memory.\n";
        log.log(LOG_WARN);
      } else {
        log << "Detected bin count " << param.num_bins << ", stddev: " << s
            << '\n';
        log.log(LOG_INFO);
      }
    }

    if (gene_matrix.n_cols < gene_matrix.n_rows) {
      log << "Your expression matrix indicates that you "
          << "have fewer genes than samples. This is rarely the case. "
          << "Please make sure genes are columns, samples are rows.\n";
      log.log(LOG_WARN);
    }

    if (param.mode == "CLR") {
      param.m = 0;
    } else if (param.mode == "ARACNE") {
      param.m = 1;
    } else if (param.mode == "RAW") {
      param.m = 2;
    }

    std::vector<std::string> targets;
    if (vm.count("targets") > 0) {
      targets = read_genes(param.targets_file);
    }

    if (vm["batch-size"].defaulted()) {
      if (vm.count("targets") > 0) {
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

    mi_full(gene_matrix, genes, targets, param);

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
