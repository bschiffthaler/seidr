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

#include <tom_fun.h>
#include <BSlogger.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#include <armadillo>
#include <boost/numeric/conversion/cast.hpp>

using boost::numeric_cast;

bool ascending(edge a, edge b) { return (a.w < b.w);}

void rank_vector(std::vector<edge>& ev)
{
  std::sort(ev.begin(), ev.end(), ascending);
  auto it = ev.begin();
  size_t pos = 0;
  double prev = it->w;
  size_t start = 0;
  double rank;
  while (it != ev.end())
  {
    it++; pos++;
    if (it == ev.end() || it->w != prev)
    {
      rank = ( numeric_cast<double>(pos) + 1 + numeric_cast<double>(start) ) / 2;
      for (size_t i = start; i < pos; i++)
      {
        edge e = ev[i];
        e.r = rank;
        ev[i] = e;
      }
      if (it != ev.end())
      {
        start = pos;
        prev = it->w;
      }
    }
  }
}


void to_rank(arma::mat& GM)
{
  GM.each_col([](arma::vec & v) {
    std::vector<edge> rank_vec;
    size_t i = 0;
    for (auto& it : v)
    {
      edge e;
      e.i = i++; e.w = (it);
      rank_vec.push_back(e);
    }
    rank_vector(rank_vec);
    for (auto& it : rank_vec)
    {
      v(it.i) = it.r;
    }
  });
}

void write_lm(arma::mat & gm, const std::string& outfile, bool abs)
{
  std::ofstream ofs(outfile, std::ios::out);
  if (abs)
  {
    gm = arma::abs(gm);
  }
  if (! ofs)
  {
    throw std::runtime_error("Could not write to file " + outfile);
  }

  for (size_t i = 1; i < gm.n_cols; i++)
  {
    for (size_t j = 0; j < i; j++)
    {
      ofs << gm(i, j);
      ofs << (j == (i - 1) ? '\n' : '\t');
    }
  }
}

double scale_free_fit(const arma::mat & gm, const uint64_t p,
                      const uint64_t nbreaks, const tom_cor_t cor)
{
  double pd = numeric_cast<double>(p);
  arma::mat gm_copy;
  if (cor == UNSIGNED)
  {
    gm_copy = arma::pow(arma::abs(gm), pd);
  }
  else if (cor == SIGNED)
  {
    gm_copy = arma::pow(((gm + 1) * 0.5), pd);
  }
  else if (cor == SIGNED_HYBRID)
  {
    gm_copy = gm;
    gm_copy.transform([](double x)->double {
      if (x < 0)
      {
        return 0;
      }
      else
      {
        return x;
      }
    });
    gm_copy = arma::pow(gm_copy, pd);
  }
  gm_copy.diag().zeros();

  // colSums
  arma::vec strengths(gm_copy.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < gm_copy.n_cols; i++)
  {
    strengths(i) = arma::accu(gm_copy.col(i));
  }

  // Vector of breaks for discretizing nodes
  double min = strengths.min();
  double max = strengths.max();
  arma::vec breaks = arma::linspace(min, max, nbreaks + 1);

  // Discretize node strengths into bins
  arma::uvec counts(nbreaks, arma::fill::zeros);
  strengths = arma::sort(strengths);
  arma::uword bin = 0;
  for (arma::uword i = 0; i < strengths.n_elem; i++)
  {
    if (bin == (counts.n_elem - 1))
    {
      counts(bin)++;
    }
    else
    {
      double s = strengths(i);
      while (s > breaks(bin + 1))
      {
        bin++;
        if (bin == (counts.n_elem - 1))
        {
          break;
        }
      }
      counts(bin)++;
    }
  }

  // Convert integer counts to doubles
  arma::vec countsd = arma::vec(nbreaks, arma::fill::zeros);
  for (arma::uword i = 0; i < counts.n_elem; i++)
  {
    countsd(i) = numeric_cast<double>(counts(i));
  }
  countsd /= arma::accu(countsd); // Convert to fractions

  // Calculate mid points of bins
  arma::vec midpoints(nbreaks, arma::fill::zeros);
  for (arma::uword i = 0; i < (breaks.n_elem - 1); i++)
  {
    midpoints(i) = ((breaks(i) + breaks(i + 1)) / 2.0);
  }

  // R^2 is the correlation of P(k) ~ k
  return arma::as_scalar(arma::abs(arma::cor(arma::log(countsd), arma::log(midpoints))));
}

arma::mat tom_similarity(arma::mat & gm, const uint64_t p,
                         const tom_cor_t cor)
{
  double pd = numeric_cast<double>(p);
  if (cor == UNSIGNED)
  {
    gm = arma::pow(arma::abs(gm), pd);
  }
  else if (cor == SIGNED)
  {
    gm = arma::pow(((gm + 1) * 0.5), pd);
  }
  else if (cor == SIGNED_HYBRID)
  {
    gm.transform([](double x)->double {
      if (x < 0)
      {
        return 0;
      }
      else
      {
        return x;
      }
    });
    gm = arma::pow(gm, pd);
  }
  gm.diag().zeros();
  arma::mat L = gm * gm;
  arma::vec K(gm.n_cols, arma::fill::zeros);
  arma::mat TOM(gm.n_rows, gm.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < gm.n_cols; i++)
  {
    K(i) = arma::accu(gm.col(i));
  }
  for (arma::uword i = 1; i < gm.n_cols; i++)
  {
    for (arma::uword j = 0; j < i; j++)
    {
      double t = (L(i, j) + gm(i, j)) / (fmin(K(i), K(j)) + 1 - gm(i, j));
      TOM(i, j) = t;
      TOM(j, i) = t;
    }
  }
  return TOM;
}

double mad(const arma::vec & v)
{
  arma::vec vp(v.n_elem, arma::fill::zeros);
  double med = arma::median(v);
  for (arma::uword i = 0; i < v.n_elem; i++)
  {
    vp(i) = (v(i) - med);
  }
  return arma::median(arma::abs(vp));
}

arma::mat bicor(const arma::mat & gm)
{
  logger log(std::cerr, "bicor");
  arma::mat gm_copy = gm;
  gm_copy.each_col([](arma::vec & col) {
    double ma = mad(col);
    double me = arma::median(col);
    arma::vec w(col.n_elem, arma::fill::zeros);
    for (arma::uword i = 0; i < col.n_elem; i++)
    {
      double xi = ((col(i) - me) / (9 * ma));
      if (fabs(xi) > 1)
      {
        w(i) = 0;
      }
      else
      {
        double ux = (1 - (xi * xi));
        w(i) = ux * ux;
      }
    }
    arma::vec denoms = col - me;
    for (arma::uword i = 0; i < denoms.n_elem; i++)
    {
      denoms(i) = denoms(i) * w(i);
    }
    double denom = sqrt(arma::accu(arma::pow(denoms, 2)));
    for (arma::uword i = 0; i < col.n_elem; i++)
    {
      col(i) =  ((col(i) - me) * w(i) / denom);
    }
  });

  arma::mat ret = gm_copy.t() * gm_copy;

  if (! ret.is_finite())
  {
    log(LOG_WARN) << "Zero valued MADs in data. Falling back to pearson "
                  << "correlation for affected columns\n";
    arma::mat cmret = arma::cor(gm);
    arma::uvec nf = arma::find_nonfinite(ret);
    ret.elem(nf) = cmret.elem(nf);
  }

  return ret;
}