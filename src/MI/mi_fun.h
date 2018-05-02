#pragma once

#include <armadillo>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <cmath>
#include <string>

class aranode {
 public:
 aranode(arma::uword a, double b) : i(a), v(b) {}
  arma::uword i;
  double v;
};

arma::vec knot_vector(size_t spline_order, size_t num_bins);
double percentile(arma::vec& data, size_t percentile);
double iqr(arma::vec& data);
double bin_width(arma::vec& data);
arma::uvec bin_count(arma::mat& gm, size_t multiplier);
arma::vec to_z(arma::vec& x);
double basis_function(size_t i, size_t p, double t,
                      arma::vec& knot_vector, size_t num_bins);
void find_weights(arma::mat& gm, arma::vec& knots, arma::mat& wm,
                  size_t spline_order, size_t num_bins, size_t i);
arma::vec hist1d(arma::vec& x, arma::vec& knots, arma::vec& weights,
                 size_t spline_order, size_t num_bins);
double log2d(double x);
double entropy1d(arma::mat& gm, arma::vec& knots, arma::mat& wm,
                 size_t spline_order, size_t num_bins, size_t i);
void hist2d(arma::vec& x, arma::vec& y, arma::vec& knots,
           arma::vec& wx, arma::vec& wy, arma::mat& hist,
            size_t spline_order, size_t num_bins);
double entropy2d(arma::mat& gm, arma::vec& knots,
                 arma::mat& wm, size_t spline_order,
                 size_t num_bins, size_t xi, size_t yi);
void mi_sub_matrix(arma::mat& gm, size_t num_bins, size_t spline_order,
                   std::vector<arma::uword> targets, std::string tempfile);
void mi_full(arma::mat& gm, size_t spline_order, size_t num_bins, size_t bs,
              std::string outfile, char m, std::string mi_file);
