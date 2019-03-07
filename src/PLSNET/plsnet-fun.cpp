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
#include <plsnet-fun.h>
#include <mpiomp.h>
#include <fs.h>
// External
#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <map>
#include <armadillo>
#include <ctime>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace fs = boost::filesystem;

std::random_device rd;
std::mt19937 gen(rd());

class seidr_mpi_plsnet : public seidr_mpi_omp {
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_predictor_sample_size(arma::uword x) {_predictor_sample_size = x;}
  void set_ensemble_size(arma::uword x) {_ensemble_size = x;}
  void set_ncomp(arma::uword x) {_ncomp = x;}
  void set_targeted(bool x) {_targeted = x;}
private:
  arma::uword _predictor_sample_size = 0;
  arma::uword _ensemble_size = 0;
  arma::uword _ncomp = 0;
  bool _targeted = 0;
};

void seidr_mpi_plsnet::entrypoint()
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  while (! _my_indices.empty())
  {
    std::vector<arma::uword> uvec;
    for (auto i : _my_indices)
    {
      uvec.push_back(i);
    }
    plsnet(_data, _genes, uvec, _tempdir, _predictor_sample_size,
           _ensemble_size, _ncomp);
    get_more_work();
  }
  log << "No more work. Waiting for other tasks to finish...\n";
  log.send(LOG_INFO);
}

void seidr_mpi_plsnet::finalize()
{
  remove(_queue_file);
  merge_files(_outfile, _tempdir,
              _targeted, _id, _genes);
}


void plsnet(const arma::mat& geneMatrix,
            const std::vector<std::string>& genes,
            const std::vector<arma::uword>& uvec,
            const std::string& tmpdir,
            const arma::uword& predictor_sample_size,
            const arma::uword& ensemble_size,
            const arma::uword& ncomp)
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  std::string tmpfile = tempfile(tmpdir);
  std::ofstream ofs(tmpfile.c_str(), std::ios::out);

  if (! ofs)
  {
    throw std::runtime_error("Could not open temp file: " + tmpfile);
  }

  // Create random generators for samples and genes
  std::uniform_int_distribution<> sample_gen(0, geneMatrix.n_rows - 1);

  arma::vec ret(geneMatrix.n_cols, arma::fill::zeros);

  #pragma omp parallel for
  for (uint64_t i = 0; i < uvec.size(); i++)
  {
    auto& target = uvec[i];
    #pragma omp critical
    {
      log << "Started gene: " << genes[target] << ".\n";
      log.send(LOG_INFO);
    }
    arma::uvec pred(geneMatrix.n_cols - 1);
    arma::uword j = 0;
    for (arma::uword i = 0; i < geneMatrix.n_cols; i++)
    {
      if (i != target)
      {
        pred(j) = i;
        j++;
      }
    }

    ret.zeros();
    std::seed_seq seeds{3, 1, 4, 1, 5, 9, 2, 6, 5};
    gen.seed(seeds);
    arma::arma_rng::set_seed(314159265);

    for (arma::uword i = 0; i < ensemble_size; i++)
    {
      // Generate vector of random samples with replacement
      arma::uvec samples(geneMatrix.n_rows);
      for (arma::uword j = 0; j < geneMatrix.n_rows; j++)
      {
        samples(j) = sample_gen(gen);
      }

      // Generate vector of random predictors without replacement
      arma::uvec pred_sub(predictor_sample_size);
      pred = arma::shuffle(pred);
      for (arma::uword j = 0; j < predictor_sample_size; j++)
      {
        pred_sub(j) = pred(j);
      }

      // Subset expression matrix
      arma::mat X = geneMatrix.submat(samples, pred_sub);
      arma::vec Y = geneMatrix.col(target);
      Y = Y(samples);
      // Get variable importance
      arma::vec vim = vip(X, Y, ncomp);
      ret(pred_sub) += vim;
    }
    #pragma omp critical
    {
      ofs << target << '\n';
      for (seidr_uword_t i = 0; i < ret.size(); i++)
      {
        ofs << ret[i] << (i == ret.size() - 1 ? '\n' : '\t');
      }
    }
  }
  ofs.close();
}

void plsnet_full(const arma::mat& GM,
                 const std::vector<std::string>& genes,
                 const seidr_plsnet_param_t& param) {
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> uvec;
  for (uint64_t i = 0; i < GM.n_cols; i++)
  {
    uvec.push_back(i);
  }

  seidr_mpi_plsnet mpi(param.bs, GM, uvec, genes, param.tempdir,
                       param.outfile);
  mpi.set_predictor_sample_size(param.predictor_sample_size);
  mpi.set_ensemble_size(param.ensemble_size);
  mpi.set_ncomp(param.ncomp);

  mpi.entrypoint();

  MPI_Barrier(MPI_COMM_WORLD); // NOLINT

  #pragma omp critical
  {
    if (mpi.rank() == 0)
    {
      while (mpi.check_logs(LOG_NAME"@" + mpi_get_host())); // NOLINT
      log << "Finalizing...\n";
      log.send(LOG_INFO);
      mpi.finalize();
    }
  }

  MPI_Finalize();
}

void plsnet_partial(const arma::mat& GM,
                    const std::vector<std::string>& genes,
                    const std::vector<std::string>& targets,
                    const seidr_plsnet_param_t& param) {

  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> positions;
  for (uint64_t i = 0; i < targets.size(); i++)
  {
    uint64_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size())
    {
      log << "Gene " << targets[i]
          << " was not found in the expression set "
          << "and will therefore not be considered."
          << " Please check that your expression set and "
          << "its column names (gene file) contain an entry for "
          << targets[i] << ".\n";
      log.log(LOG_WARN);
    }
    else
    {
      positions.push_back(pos);
    }
  }

  seidr_mpi_plsnet mpi(param.bs, GM, positions, genes, param.tempdir,
                       param.outfile);
  mpi.set_predictor_sample_size(param.predictor_sample_size);
  mpi.set_ensemble_size(param.ensemble_size);
  mpi.set_ncomp(param.ncomp);
  mpi.set_targeted(true);

  mpi.entrypoint();

  MPI_Barrier(MPI_COMM_WORLD); // NOLINT

  #pragma omp critical
  {
    if (mpi.rank() == 0)
    {
      while (mpi.check_logs(LOG_NAME"@" + mpi_get_host())); // NOLINT
      log << "Finalizing...\n";
      log.send(LOG_INFO);
      mpi.finalize();
    }
  }

  MPI_Finalize();
}

arma::vec vip(arma::mat& X, arma::vec& Y, arma::uword ncomp)
{
  plsreg_t plsr = plsreg(X, Y, ncomp);
  arma::mat wsq = arma::pow(plsr.w, 2);
  wsq = wsq.t();
  double xnr = X.n_cols;
  arma::mat xret = xnr * plsr.p.row(1) * wsq / arma::accu(plsr.p.row(1));
  xret = xret.t();
  return xret.col(0);
}

plsreg_t plsreg(arma::mat& X, arma::vec& Y, arma::uword ncomp)
{
  plsreg_t ret;

  arma::mat X0 = X;
  X0.each_col([](arma::vec & v) { v -= arma::mean(v); });

  arma::vec Y0 = Y - arma::mean(Y);

  simpls_t spls = simpls(X0, Y0, ncomp);

  arma::mat pctvar(2, ncomp, arma::fill::zeros);

  pctvar.row(0) = arma::sum(arma::pow(spls.x, 2), 0) /
                  arma::accu(arma::sum(arma::pow(X0, 2), 0));
  pctvar.row(1) = arma::sum(arma::pow(spls.y, 2), 0) /
                  arma::accu(arma::sum(arma::pow(Y0, 2), 0));

  ret.p = pctvar;
  ret.w = spls.w;

  return ret;
}

simpls_t simpls(arma::mat& X, arma::vec& Y, arma::uword ncomp)
{
  simpls_t ret;
  arma::mat Cov = X.t() * Y;

  arma::mat Xloadings(X.n_cols, ncomp, arma::fill::zeros);
  arma::rowvec Yloadings(ncomp, arma::fill::zeros);
  arma::mat Weights = Xloadings;
  arma::mat W = Xloadings;
  //cov.print();
  for (arma::uword i = 0; i < ncomp; i++)
  {
    arma::mat U, V;
    arma::vec s;
    svd_econ(U, s, V, Cov);
    U = U.col(0);
    s = s(0);
    V = V.col(0);
    arma::vec ti = X * U;
    double normti = norm(ti);
    ti /= normti;
    Xloadings.col(i) = X.t() * ti;
    double qi = s(0) * V(0, 0) / normti;
    Yloadings(i) = qi;
    U /= normti;
    Weights.col(i) = U;
    arma::vec vi = Xloadings.col(i);
    for (arma::uword r = 0; r < 2; r++)
    {
      for (arma::uword j  = 0; j < i; j++)
      {
        arma::vec vj = W.col(j);
        arma::vec tmp = (vi.t() * vj);
        vi -= vj * tmp(0);
      };
    }
    vi /= norm(vi);
    W.col(i) = vi;
    Cov = Cov - vi * (vi.t() * Cov);
    arma::mat Vi = W.cols(0, i);
    Cov = Cov - Vi * (Vi.t() * Cov);
  }
  ret.x = Xloadings;
  ret.y = Yloadings;
  ret.w = Weights;
  return ret;
}
