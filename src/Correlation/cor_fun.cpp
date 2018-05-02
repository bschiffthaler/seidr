#include <armadillo>
#include <vector>
#include <algorithm>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <cor_fun.h>

using boost::lexical_cast;

bool ascending(edge a, edge b){ return (a.w < b.w);}

void rank_vector(std::vector<edge>& ev)
{
  std::sort(ev.begin(), ev.end(), ascending);
  auto it = ev.begin();
  size_t pos = 0;
  double prev = it->w;
  size_t start = 0;
  double rank;
  while(it != ev.end())
    {
      it++; pos++;
      if(it->w != prev || it == ev.end())
        {
          rank = ( lexical_cast<double>(pos) + 1 + lexical_cast<double>(start) ) / 2;
          for(size_t i = start; i < pos; i++)
            {
              edge e = ev[i];
              e.r = rank;
              ev[i] = e;
            }
          if(it != ev.end())
            {
              start = pos;
              prev = it->w;
            }
        }
    }
}


void to_rank(arma::mat& GM)
{
  GM.each_col([](arma::vec& v){
      std::vector<edge> rank_vec;
      size_t i = 0;
      for(auto it = v.begin(); it != v.end(); it++)
        {
          edge e;
          e.i = i++; e.w = (*it);
          rank_vec.push_back(e);
        }
      rank_vector(rank_vec);
      for(auto it = rank_vec.begin(); it != rank_vec.end(); it++)
        {
          v(it->i) = it->r;
        }
    });
}

void write_lm(arma::mat gm, std::string outfile, bool abs)
{
  std::ofstream ofs(outfile, std::ios::out);
  if(abs)
    gm = arma::abs(gm);
  if(! ofs)
    throw std::runtime_error("Could not write to file " + outfile);
  
  for(size_t i = 1; i < gm.n_cols; i++)
    {
      for(size_t j = 0; j < i; j++)
        {
          ofs << gm(i,j);
          ofs << (j == (i - 1) ? '\n' : '\t');
        }
    }
  
}
