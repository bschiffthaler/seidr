#pragma once

#include <common.h>
#include <armadillo>
#include <vector>
#include <string>

void genie3(arma::mat gm, std::vector<std::string> genes,
            std::vector<arma::uword> pred, std::string outfile,
            int thread_id,
            seidr_uword_t ntree, seidr_uword_t mtry,
            size_t min_node_size, double alpha, double minprop);
void genie3_full(arma::mat gm, std::vector<std::string> genes,
                 arma::uword bs, std::string outfile,
                 seidr_uword_t ntree, seidr_uword_t mtry,
                 size_t min_node_size, double alpha, double minprop);
void genie3_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                    std::vector<std::string> targets, std::string outfile,
                    seidr_uword_t nstree, seidr_uword_t mtry,
                    size_t min_node_size, double alpha, double minprop);
