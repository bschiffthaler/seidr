/* * Seidr - Create and operate on gene crowd networks
 * Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
 *
 * This file is part of Seidr.
 *
 * Seidr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Seidr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include <armadillo>
// covariance between two vectors
arma::mat
mcov(const arma::colvec& x, const arma::colvec& y);
// calculate mutual information
double
cmi(const arma::colvec& x, const arma::colvec& y);
// calculate MI for a whole set of predictors
arma::vec
fullMI(const arma::mat& X, const arma::colvec& Y);
// cumulative distribution function(x)
double
phi(const double& x, const double& mu, const double& sigma);

arma::vec
normcdf(const arma::vec& x, const double& mu, const double& sigma);
