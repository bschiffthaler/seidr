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

#include <cor_fun.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <armadillo>
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;

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
      rank = ( lexical_cast<double>(pos) + 1 + lexical_cast<double>(start) ) / 2;
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

void write_lm(arma::mat& gm, const std::string& outfile, bool abs)
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
