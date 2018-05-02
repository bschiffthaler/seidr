#pragma once

#include <vector>
#include <armadillo>
#include <string>

std::vector<arma::uvec> shuffle(unsigned int);
void tiglm(arma::mat, std::vector<std::string>,
           std::vector<arma::uword>, std::string, int,
           seidr_uword_t boot, double fmin,
           seidr_uword_t nsteps);
void tiglm_full(arma::mat, std::vector<std::string>, size_t, std::string,
                seidr_uword_t boot, double fmin,
                seidr_uword_t nsteps);
void tiglm_partial(arma::mat, std::vector<std::string>, size_t, 
		   std::vector<std::string>, std::string,
                   seidr_uword_t boot, double fmin,
                   seidr_uword_t nsteps);
