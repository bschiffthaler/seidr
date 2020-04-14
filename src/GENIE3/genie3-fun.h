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

#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif
#include <armadillo>
#include <common.h>
#include <string>
#include <vector>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#define GENIE3_FULL 0
#define GENIE3_PARTIAL 1

#define LOG_NAME "genie3"

class seidr_mpi_genie3;

struct seidr_genie3_param_t
{
  friend class boost::serialization::access;
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& BOOST_SERIALIZATION_NVP(infile);
    ar& BOOST_SERIALIZATION_NVP(gene_file);
    ar& BOOST_SERIALIZATION_NVP(targets_file);
    ar& BOOST_SERIALIZATION_NVP(do_scale);
    ar& BOOST_SERIALIZATION_NVP(force);
    ar& BOOST_SERIALIZATION_NVP(row_delim);
    ar& BOOST_SERIALIZATION_NVP(field_delim);
    ar& BOOST_SERIALIZATION_NVP(bs);
    ar& BOOST_SERIALIZATION_NVP(mode);
    ar& BOOST_SERIALIZATION_NVP(outfile);
    ar& BOOST_SERIALIZATION_NVP(tempdir);
    ar& BOOST_SERIALIZATION_NVP(verbosity);
    ar& BOOST_SERIALIZATION_NVP(nthreads);
    ar& BOOST_SERIALIZATION_NVP(ntree);
    ar& BOOST_SERIALIZATION_NVP(mtry);
    ar& BOOST_SERIALIZATION_NVP(min_node_size);
    ar& BOOST_SERIALIZATION_NVP(minprop);
    ar& BOOST_SERIALIZATION_NVP(alpha);
  }
  std::string infile;
  std::string gene_file;
  std::string targets_file;
  bool do_scale = false;
  bool force = false;
  char row_delim = '\n';
  char field_delim = '\t';
  uint64_t bs;
  uint64_t mode = GENIE3_FULL;
  std::string outfile;
  uint64_t ntree;
  uint64_t mtry;
  uint64_t min_node_size;
  double alpha;
  double minprop;
  std::string tempdir;
  unsigned verbosity;
  bool resuming;
  int nthreads;
  std::vector<uint64_t> good_idx;
  std::string cmd_file;
};

void
genie3(const arma::mat& gm,
       const std::vector<std::string>& genes,
       const std::vector<arma::uword>& pred,
       const std::string& tmpdir,
       const uint64_t& ntree,
       const uint64_t& mtry,
       const uint64_t& min_node_size,
       const double& alpha,
       const double& minprop,
       seidr_mpi_genie3* self);

void
genie3_full(const arma::mat& gm,
            const std::vector<std::string>& genes,
            const seidr_genie3_param_t& param);

void
genie3_partial(const arma::mat& gm,
               const std::vector<std::string>& genes,
               const std::vector<std::string>& targets,
               const seidr_genie3_param_t& param);
