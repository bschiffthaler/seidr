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
#include <BSlogger.hpp>
#include <IP_LPT2.h>
#include <common.h>
// External
#include <armadillo>
#include <boost/numeric/conversion/cast.hpp>
#include <coin/ClpInterior.hpp>
#include <coin/ClpSimplex.hpp>
#include <coin/ClpSolve.hpp>
#include <coin/CoinPackedMatrix.hpp>
#include <omp.h>
#include <vector>

// Solve linear programming optimisation
arma::mat
lp_ipt(const arma::mat& X, const arma::rowvec& Y, const std::string& al)
{
  uint64_t ng = X.n_rows; // number of genes
  uint64_t ns = X.n_cols; // number of samples
  uint64_t lps =
    2 * ns + 2 * X.n_elem + 1; // size of linear problem matrix vector
  arma::mat xf = arma::join_horiz(X.t(), X.t() * -1); // right side matrix
  uint64_t pc = 2 * ns + xf.n_cols; // columns of problem matrix

  // Row index = ia, Col index ij; value ir
  // std::cerr << lps << std::endl;

  std::vector<int> ia(lps, 0);
  std::vector<int> ja(lps, 0);
  std::vector<double> ar(lps, 0);

  std::vector<double> obj(pc, 0);
  std::vector<double> rowlb(Y.n_elem, 0);
  std::vector<double> rowub(Y.n_elem, 0);

  for (uint64_t i = 0; i < pc; i++) {
    obj[i] = 1;
  }

  for (uint64_t i = 0; i < Y.n_elem; i++) {
    rowlb[i] = Y(i);
    rowub[i] = Y(i);
  }

  for (uint64_t i = 0; i < ns; i++) {
    ia[i] = i;
    ja[i] = i;
    ar[i] = 1;

    ia[i + ns] = i;
    ja[i + ns] = i + ns;
    ar[i + ns] = -1;
  }

  uint64_t vo = 2 * ns; // value offset in matrix vectors
  for (uint64_t i = 0; i < xf.n_rows; i++) {
    for (uint64_t j = 0; j < xf.n_cols; j++) {
      ia[vo] = i;
      ja[vo] = j + 2 * ns;
      ar[vo] = xf(i, j);
      vo++;
    }
  }

  CoinPackedMatrix problem_mat(
    true, ia.data(), ja.data(), ar.data(), boost::numeric_cast<int>(lps));

  double const* solution;
  if (al == "interior-point") {
    auto model = ClpInterior();
    model.setLogLevel(0); // Be quiet
    model.loadProblem(
      problem_mat, NULL, NULL, obj.data(), rowlb.data(), rowub.data());
    model.setOptimizationDirection(1); // Minimize
    model.setPrimalTolerance(1e-7);    // Same as GLPK
    model.primalDual();
    solution = model.primalColumnSolution();
  } else {
    auto model = ClpSimplex();
    model.setLogLevel(0); // Be quiet
    model.loadProblem(
      problem_mat, NULL, NULL, obj.data(), rowlb.data(), rowub.data());
    model.setOptimizationDirection(1); // Minimize
    model.setPrimalTolerance(1e-7);    // Same as GLPK
    model.initialSolve();
    solution = model.primalColumnSolution();
  }

  arma::vec res(xf.n_cols + 2 * ns);
  for (uint64_t i = 0; i < xf.n_cols + 2 * ns; i++) {
    res(i) = solution[i];
  }

  arma::mat fin(1, ng);
  for (uint64_t i = 0; i < ng; i++) {
    fin(0, i) = res(2 * ns + i);
    fin(0, i) -= res(2 * ns + ng + i);
  }

  return fin;
}
