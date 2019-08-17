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
#include <narromi_fun.h>
#include <common.h>
#include <IP_LPT.h>
#include <stats_fun.h>
#include <mpiomp.h>
#include <fs.h>
// External
#include <mpi.h>
#include <armadillo>
#include <vector>
#include <fstream>
#include <sys/types.h>
#include <sys/wait.h>
#include <ctime>
#include <boost/filesystem.hpp>
#include <map>

namespace fs = boost::filesystem;

NaResult::NaResult(arma::vec sig, arma::vec nv,
                   arma::uvec re, arma::uword i) {
  s = sig; n = nv; r = re; I = i;
}
arma::vec NaResult::sig() { return s; }
arma::vec NaResult::nv() { return n; }
arma::uvec NaResult::re() { return r; }
arma::uword NaResult::i() { return I; }


class seidr_mpi_narromi : public seidr_mpi_omp {
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_algorithm(std::string a) {_algorithm = a;}
  void set_alpha(double a) {_alpha = a;}
  void set_t(double t) {_t = t;}
  void set_targeted(bool x) {_targeted = x;}
private:
  std::string _algorithm = "simplex";
  double _alpha = 0;
  double _t = 0;
  bool _targeted = false;
};

void seidr_mpi_narromi::entrypoint()
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  while (! _my_indices.empty())
  {
    std::vector<arma::uword> uvec;
    for (auto i : _my_indices)
    {
      uvec.push_back(i);
    }
    narromi_thread(_data, _algorithm, _alpha, _t, uvec, _genes, _tempdir);
    get_more_work();
  }
  #pragma omp critical
  {
    log << "No more work. Waiting for other tasks to finish...\n";
    log.send(LOG_INFO);
  }
}

void seidr_mpi_narromi::finalize()
{
  merge_files(_outfile, _tempdir,
              _targeted, _id, _genes);
}

/* Linear programming optimisation loop.
 */
std::pair<arma::vec, arma::vec>
reoptim(const arma::mat& X,
        const arma::rowvec& Y,
        const std::string& al,
        const double& alpha) {

  arma::mat J = lp_ipt(X, Y, al);
  arma::vec Js(J.n_cols);
  arma::vec Jv(J.n_cols);

  Js.zeros(); Jv.zeros();

  arma::uvec index   = arma::find(arma::abs(J) >= alpha);
  arma::uvec index_c = arma::find(arma::abs(J) < alpha);

  Js(index) = J(index);
  Jv(index_c) = J(index_c);

  while (index_c.n_elem > 0) {
    arma::mat J1 = lp_ipt(X.rows(index), Y, al);
    arma::uvec index1   = arma::find(arma::abs(J1) >= alpha);
    arma::uvec index1_c = arma::find(arma::abs(J1) < alpha);
    index_c     = index(index1_c);
    index       = index(index1);
    Js(index)   = J1(index1);
    Js(index_c).zeros();
    Jv(index_c) = J1(index1_c);
  }

  Jv = Js + Jv;

  return std::pair<arma::vec, arma::vec> (Js, Jv);
}

/* Significance calculation function */

arma::vec get_sig(arma::vec z) {
  z = (z - arma::min(z)) / (arma::max(z) - arma::min(z));
  z = 0.5 * log2((1 + z) / (1 - z));
  // Set max values to second highest values
  arma::uvec m1 = arma::find(z == arma::max(z));
  arma::uvec m2 = arma::find(z != arma::max(z));
  z(m1).fill(arma::max(z(m2)));
  double mu = mean(z);
  double sigma = var(z);
  arma::vec sig = 1 - normcdf(z, mu, sigma);
  return sig;
}

// Run narromi algorithm for one gene against a set
// of genes
NaResult narromi(const arma::mat& GM,
                 const std::string& al,
                 const double& alpha,
                 const double& t,
                 const arma::uword& i) {

  arma::uvec response = get_i(i, GM.n_cols);
  arma::mat X = GM.cols(response);
  arma::vec Y = GM.col(i);
  arma::vec net(X.n_cols, arma::fill::zeros);
  arma::vec fmi       = fullMI(X, Y);
  arma::vec net_value = fmi;
  arma::uvec index    = find(arma::abs(fmi) >= alpha);
  arma::mat X1        = X.cols(index);

  std::pair<arma::vec, arma::vec> J = reoptim(X1.t(), Y.t(), al, alpha);

  for (size_t i = 0; i < index.n_elem; i++) {
    net(index)        = J.first;
    net_value(index)  = J.second;
  }

  arma::vec net_value1 = net_value;
  net_value = arma::sign(net_value) % (arma::abs(net_value) * t + fmi * (1 - t));
  arma::vec z = arma::abs(net_value1);
  arma::vec sig1 = get_sig(z);

  z = (arma::abs(fmi));
  arma::vec sig2 = get_sig(z);

  arma::vec sig = arma::sqrt(arma::pow(sig1, 2) + arma::pow(sig2, 2));

  return NaResult(sig, net_value, response, i);

}

/* Run narromi algorithm on a full expression set in an all vs.
all style comparison.
 */
void full_narromi(const arma::mat& GM,
                  const std::vector<std::string>& genes,
                  const seidr_narromi_param_t& param) {
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> uvec;
  for (uint64_t i = 0; i < GM.n_cols; i++)
  {
    uvec.push_back(i);
  }

  seidr_mpi_narromi mpi(param.bs, GM, uvec, genes, param.tempdir,
                        param.outfile);
  mpi.set_alpha(param.alpha);
  mpi.set_t(param.t);
  mpi.set_algorithm(param.al);
  mpi.entrypoint();

  MPI_Barrier(MPI_COMM_WORLD); // NOLINT
  mpi.remove_queue_file();
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

/* Run narromi algorithm on a full expression set on a set of
genes of interest.
*/

void partial_narromi(const arma::mat& GM,
                     const std::vector<std::string>& genes,
                     const std::vector<std::string>& targets,
                     const seidr_narromi_param_t& param) {

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  std::vector<uint64_t> positions;
  for (uint64_t i = 0; i < targets.size(); i++)
  {
    uint64_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size())
    {
      log << "\nWarning: Gene " << targets[i]
          << " was not found in the expression set "
          << "and will therefore not be considered. "
          << "Please check that your expression set and "
          << "its column (gene) names contain an entry for " << targets[i] << ".\n";
      log.log(LOG_WARN);
    }
    else
    {
      positions.push_back(pos);
    }
  }

  seidr_mpi_narromi mpi(param.bs, GM, positions, genes, param.tempdir,
                        param.outfile);
  mpi.set_alpha(param.alpha);
  mpi.set_t(param.t);
  mpi.set_algorithm(param.al);
  mpi.set_targeted(true);
  mpi.entrypoint();

  MPI_Barrier(MPI_COMM_WORLD); // NOLINT
  mpi.remove_queue_file();
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

void narromi_thread(const arma::mat& gene_matrix,
                    const std::string& al,
                    const double& alpha,
                    const double& t,
                    const std::vector<arma::uword>& ind,
                    const std::vector<std::string>& genes,
                    const std::string tmpdir) {
  std::clock_t start;
  size_t tot = gene_matrix.n_cols - 1;
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  std::string tmpfile = tempfile(tmpdir);
  std::ofstream ofs(tmpfile.c_str());

  if (! ofs)
  {
    throw std::runtime_error("Could not open temp file: " + tmpfile);
  }

  #pragma omp parallel for
  for (size_t x = 0; x < ind.size(); x++) 
  {
    start = std::clock();
    arma::uword i = ind[x];
    #pragma omp critical
    {
      log << "Started gene " << genes[i] << '\n';
      log.send(LOG_INFO);
    }
    NaResult nr = narromi(gene_matrix, al, alpha, t, i);
    arma::uvec re = nr.re();
    arma::vec si = nr.sig();
    arma::vec nv = nr.nv();

    #pragma omp critical
    {
      ofs << i << '\n';
      for (arma::uword j = 0; j < tot; j++) {
        if (j == i && j < (tot - 1))
        {
          ofs << 0 << '\t' << nv(j) << '\t';
        }
        else if (i == j && j == (tot - 1))
        {
          ofs << 0 << '\t' << nv(j) << '\n';
        }
        else if (i > j && j == (tot - 1))
        {
          ofs << nv(j) << '\t' << 0 << '\n';
        }
        else
        {
          ofs << nv(j) <<
              (j == (tot - 1) ? '\n' : '\t');
        }
      }
    }
  }
  ofs.close();
}
