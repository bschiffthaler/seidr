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
#include <glmnet2.h>
#include <tiglm.h>
#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif
// External
#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>

namespace fs = boost::filesystem;

class seidr_mpi_tigress : public seidr_mpi_omp
{
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_nboot(seidr_uword_t x) { _nboot = x; }
  void set_fmin(double x) { _fmin = x; }
  void set_nsteps(seidr_uword_t x) { _nsteps = x; }
  seidr_uword_t get_nboot() const { return _nboot; }
  double get_fmin() const { return _fmin; }
  seidr_uword_t set_nsteps() const { return _nsteps; }
  void set_targeted(bool x) { _targeted = x; }
  void set_allow_low_var(bool x) { _allow_low_var = x; }

private:
  seidr_uword_t _nboot = 0;
  seidr_uword_t _nsteps = 0;
  double _fmin = 0;
  bool _targeted = false;
  bool _allow_low_var = false;
};

void
seidr_mpi_tigress::entrypoint()
{
  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());
  while (!_my_indices.empty()) {
    std::vector<arma::uword> uvec;
    for (auto i : _my_indices) {
      uvec.push_back(i);
    }
    while (check_logs(LOG_NAME "@" + mpi_get_host())) {}
    tiglm(_data,
          _genes,
          uvec,
          _tempdir,
          _nboot,
          _fmin,
          _nsteps,
          _allow_low_var,
          this);
    get_more_work();
  }
}

void
seidr_mpi_tigress::finalize()
{
  merge_files(_outfile, _tempdir, _targeted, _id, _genes);
}

/*
In order to do stability selection we need to split the samples
in half at random. This function takes a sample size and outputs
a std::vector containing two arma::uvec vectors which hold the
randomized indices of the matrix rows.
*/
std::vector<arma::uvec>
shuffle(unsigned int n)
{
  arma::uvec samples(n);
  unsigned int nA = 0;
  unsigned int nB = 0;

  // Determine if sample size is even or uneven to get the right
  // size for the uvec objects
  if (n % 2 == 0) {
    nA = n / 2;
    nB = n / 2;
  } else {
    unsigned int nT = n + 1;
    nA = nT / 2;
    nB = nT / 2 - 1;
  }
  arma::uvec resA(nA);
  arma::uvec resB(nB);
  std::vector<arma::uvec> final;

  for (unsigned int i = 0; i < n; i++) {
    samples[i] = i;
  }

  arma::uvec samplesS = arma::shuffle(samples);

  for (unsigned int i = 0; i < n; i++) {
    if (i < nA) {
      resA[i] = samplesS[i];
    }
    if (i >= nA) {
      resB[i - nA] = samplesS[i];
    }
  }

  final.push_back(resA);
  final.push_back(resB);
  return final;
}

void
tiglm(const arma::mat& geneMatrix,
      const std::vector<std::string>& genes,
      const std::vector<arma::uword>& pred,
      const std::string& tmpdir,
      const seidr_uword_t& boot,
      const double& fmin,
      const seidr_uword_t& nsteps,
      const bool& allow_low_var,
      seidr_mpi_tigress* self)
{
  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());

  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> distribution(TIGRESS_SAMPLE_MIN, TIGRESS_SAMPLE_MAX);


  size_t ng = geneMatrix.n_cols - 1;
  size_t ns = geneMatrix.n_rows;
  double ns_d = geneMatrix.n_rows;

  size_t n_good = ceil(sqrt(ns_d / 2));
  if (allow_low_var) {
    n_good = 1;
  }

  std::string tmpfile = tempfile(tmpdir);
  std::ofstream ofs(tmpfile.c_str(), std::ios::out);

  if (!ofs) {
    throw std::runtime_error("Could not open temp file: " + tmpfile);
  }

#pragma omp parallel for schedule(dynamic) // NOLINT
  for (uint64_t i = 0; i < pred.size(); i++) { // NOLINT
    arma::uvec indices = get_i(pred[i], ng + 1);
    arma::mat cum = arma::zeros<arma::mat>(ng, nsteps);

    arma::mat response = geneMatrix.cols(indices);
    arma::vec predictor = geneMatrix.col(pred[i]);

    size_t iter = 0;
    size_t nboot = 0;
    /* This int is incremented if the random shuffle produced a vector
       with all 0s. If that happens too often, we skip the gene and
       return no interations, as the other case could result in an
       endless loop or a non-converging lasso. */

    arma::vec randV(ng, 1);

    std::seed_seq sseq SEIDR_DEFAULT_SSEQ;
    generator.seed(sseq);
    arma::arma_rng::set_seed(SEIDR_DEFAULT_SEED);

    while (iter < boot) {

#ifdef DEBUG
#pragma omp critical
      {
        log << "Iter: " << iter << '\n';
        log.send(LOG_DEBUG);
        while (self->check_logs(LOG_NAME "@" + mpi_get_host())) {}
      }
#endif

      // Shuffle the samples
      std::vector<arma::uvec> split = shuffle(ns);
      arma::vec yA = predictor(split[0]);
      arma::vec yB = predictor(split[1]);

      arma::uvec syA = arma::find_unique(yA);
      arma::uvec syB = arma::find_unique(yB);

      // If any of the vectors have less than two usable values, try again
      if (syA.n_elem < n_good || syB.n_elem < n_good) {
#ifdef DEBUG
#pragma omp critical
        {
          log << "Have to shuffle: A=" << syA.n_elem << " B=" << syB.n_elem
              << '\n';
          log.send(LOG_DEBUG);
          while (self->check_logs(LOG_NAME "@" + mpi_get_host())) {}
        }
#endif
        nboot++;
        if (nboot > 1000) { // NOLINT
#pragma omp critical
          {
            log << genes[pred[i]]
                << " encountered too many reshuffles... skipping.\n";
            log.send(LOG_WARN);
          }
          cum.fill(0);
          iter = boot;
        }
        continue;
      }

      randV.imbue([&]() { return distribution(generator); });
      arma::mat rw = response.t();
      for (unsigned int y = 0; y < ns; y++) {
        rw.col(y) % randV;
      }
      rw = rw.t();

      arma::mat xA = rw.rows(split[0]);
      arma::mat xB = rw.rows(split[1]);

      glm pathA = glmnet(xA, yA, nsteps, fmin);
      glm pathB = glmnet(xB, yB, nsteps, fmin);

      if (pathA.success && pathB.success) {
        cum += arma::sign(arma::abs(pathA.beta));
        cum += arma::sign(arma::abs(pathB.beta));
        iter++;
      }
    }

    auto d_nsteps = boost::numeric_cast<double>(nsteps);

    cum = arma::cumsum(cum.t()) / (boot * 2) / d_nsteps;

#pragma omp critical
    { // NOLINT
      self->increment_pbar();
      arma::uword xi = pred[i];
      ofs << xi << '\n';
      arma::rowvec fin = cum.row(nsteps - 1);
      for (arma::uword s = 0; s < indices.n_elem; s++) {
        if (s == xi && s < (indices.n_elem - 1)) {
          ofs << 0 << '\t' << fin(s) << '\t';
        } else if (xi == s && s == (indices.n_elem - 1)) {
          ofs << 0 << '\t' << fin(s) << '\n';
        } else if (xi > s && s == (indices.n_elem - 1)) {
          ofs << fin(s) << '\t' << 0 << '\n';
        } else {
          ofs << fin(s) << (s == (indices.n_elem - 1) ? '\n' : '\t');
        }
      }
    }
  }
}

void
tiglm_full(const arma::mat& GM,
           const std::vector<std::string>& genes,
           const seidr_tigress_param_t& param)
{
  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> uvec;
  for (uint64_t i = 0; i < GM.n_cols; i++) {
    if (param.resuming) {
      if (!in_sorted_range<uint64_t>(i, param.good_idx)) {
        uvec.push_back(i);
      }
    } else {
      uvec.push_back(i);
    }
  }

  seidr_mpi_tigress mpi(
    param.bs, GM, uvec, genes, param.tempdir, param.outfile);
  mpi.set_nboot(param.boot);
  mpi.set_fmin(param.fmin);
  mpi.set_nsteps(param.nsteps);
  mpi.set_allow_low_var(param.allow_low_var);

  mpi.entrypoint();

  SEIDR_MPI_BARRIER(); // NOLINT
  mpi.remove_queue_file();
#pragma omp critical
  {
    if (mpi.rank() == 0) {
      while (mpi.check_logs(LOG_NAME "@" + mpi_get_host())) {}
      mpi.finalize_pbar();
      log << "Finalizing...\n";
      log.send(LOG_INFO);
      mpi.finalize();
      while (mpi.check_logs(LOG_NAME "@" + mpi_get_host())) {}
    }
  }

  SEIDR_MPI_FINALIZE();
}

void
tiglm_partial(const arma::mat& GM,
              const std::vector<std::string>& genes,
              const std::vector<std::string>& targets,
              const seidr_tigress_param_t& param)
{

  seidr_mpi_logger log(LOG_NAME "@" + mpi_get_host());

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> positions;
  for (const auto& target : targets) {
    uint64_t pos = find(genes.begin(), genes.end(), target) - genes.begin();
    if (pos >= genes.size()) {
      log << "Gene " << target << " was not found in the expression set "
          << "and will therefore not be considered."
          << " Please check that your expression set and "
          << "its column names (gene file) contain an entry for " << target
          << ".\n";
      log.log(LOG_WARN);
    } else {
      if (param.resuming) {
        if (!in_sorted_range<uint64_t>(pos, param.good_idx)) {
          positions.push_back(pos);
        }
      } else {
        positions.push_back(pos);
      }
    }
  }

  seidr_mpi_tigress mpi(
    param.bs, GM, positions, genes, param.tempdir, param.outfile);
  mpi.set_nboot(param.boot);
  mpi.set_fmin(param.fmin);
  mpi.set_nsteps(param.nsteps);
  mpi.set_targeted(true);
  mpi.set_allow_low_var(param.allow_low_var);

  mpi.entrypoint();

  SEIDR_MPI_BARRIER(); // NOLINT
  mpi.remove_queue_file();
#pragma omp critical
  {
    if (mpi.rank() == 0) {
      while (mpi.check_logs(LOG_NAME "@" + mpi_get_host())) {}
      mpi.finalize_pbar();
      log << "Finalizing...\n";
      log.send(LOG_INFO);
      mpi.finalize();
      while (mpi.check_logs(LOG_NAME "@" + mpi_get_host())) {}
    }
  }

  SEIDR_MPI_FINALIZE();
}
