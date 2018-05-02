#pragma once
#include <armadillo>
#include <vector>
#include <string>

class NaResult {
  arma::vec s;
  arma::vec n;
  arma::uvec r;
  arma::uword I;
 public:
  NaResult(arma::vec, arma::vec,
           arma::uvec, arma::uword);
  arma::vec sig();
  arma::vec nv();
  arma::uvec re();
  arma::uword i();
};

std::pair<arma::vec, arma::vec>
  reoptim(arma::mat X, arma::rowvec Y,
	  std::string al, double alpha);

NaResult  narromi(arma::mat GM,
	      std::string al, double alpha, double t,
	      arma::uword i);

arma::vec get_sig(arma::vec z);

void full_narromi(arma::mat gene_matrix,
		  std::string al, double alpha, double t,
		  std::vector<std::string> genes, size_t nt,
                  std::string outfile);

void partial_narromi(arma::mat gene_matrix,
		     std::string al, double alpha, double t,
		     std::vector<std::string> genes, size_t nt, 
		     std::vector<std::string> targets,
                     std::string outfile);

void narromi_thread(arma::mat gene_matrix,
                    std::string al, double alpha, double t,
                    std::vector<arma::uword> ind,
                    std::vector<std::string> genes,
                    std::string tempfile, int id);

