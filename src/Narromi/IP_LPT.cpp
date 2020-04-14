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
#ifdef NARROMI_USE_CLP
#include <IP_LPT2.h>
#else
#include <IP_LPT.h>
#endif
#include <BSlogger.hpp>
#include <common.h>
// External
#include <armadillo>
#include <glpk.h>
#include <omp.h>

// Solve linear programming optimisation
arma::mat
lp_ipt(const arma::mat& X, const arma::rowvec& Y, const std::string& al)
{
  /* declare variables
   */
  glp_smcp sparm;
  glp_iptcp iparm;

  if (al == "simplex") {
    glp_init_smcp(&sparm);
    sparm.meth = GLP_DUALP;
    sparm.msg_lev = GLP_MSG_ERR;
  } else if (al == "interior-point") {
    glp_init_iptcp(&iparm);
    iparm.msg_lev = GLP_MSG_ERR;
  }

  size_t ng = X.n_rows;                   // number of genes
  size_t ns = X.n_cols;                   // number of samples
  size_t lps = 2 * ns + 2 * X.n_elem + 1; // size of linear problem matrix
                                          // vector
  arma::mat xf = arma::join_horiz(X.t(), X.t() * -1); // right side matrix
  size_t pc = 2 * ns + xf.n_cols; // columns of problem matrix

  // Row index = ia, Col index ij; value ir
  // std::cerr << lps << std::endl;
  int* ia = (int*)malloc(lps * sizeof(int));
  int* ja = (int*)malloc(lps * sizeof(int));
  double* ar = (double*)malloc(lps * sizeof(double));

  glp_prob* lp;
  lp = glp_create_prob();
  glp_set_obj_dir(lp, GLP_MIN);

  glp_add_rows(lp, Y.n_elem);
  for (size_t i = 0; i < Y.n_elem; i++) {
    glp_set_row_bnds(lp, i + 1, GLP_FX, Y(i), 0.0);
  }

  glp_add_cols(lp, pc);
  for (size_t i = 1; i <= pc; i++) {
    glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
    glp_set_obj_coef(lp, i, 1.0);
  }
  for (size_t i = 1; i <= ns; i++) {
    ia[i] = i;
    ja[i] = i;
    ar[i] = 1;

    ia[i + ns] = i;
    ja[i + ns] = i + ns;
    ar[i + ns] = -1;
  }

  size_t vo = 2 * ns + 1; // value offset in matrix vectors
  for (size_t i = 0; i < xf.n_rows; i++) {
    for (size_t j = 0; j < xf.n_cols; j++) {
      ia[vo] = i + 1;
      ja[vo] = j + 2 * ns + 1;
      ar[vo] = xf(i, j);
      vo++;
    }
  }

  glp_load_matrix(lp, lps - 1, ia, ja, ar);

  if (al == "interior-point")
    glp_interior(lp, &iparm);
  if (al == "simplex")
    // glp_exact(lp, &sparm);
    glp_simplex(lp, &sparm);

  arma::vec res(xf.n_cols + 2 * ns);
  for (size_t i = 0; i < xf.n_cols + 2 * ns; i++) {
    if (al == "interior-point")
      res(i) = glp_ipt_col_prim(lp, i + 1);
    if (al == "simplex")
      res(i) = glp_get_col_prim(lp, i + 1);
  }

  arma::mat fin(1, ng);
  for (size_t i = 0; i < ng; i++) {
    fin(0, i) = res(2 * ns + i);
    fin(0, i) -= res(2 * ns + ng + i);
  }

  /* housekeeping */
  glp_delete_prob(lp);
  glp_free_env();
  free(ia);
  free(ja);
  free(ar);
  return fin;
}
