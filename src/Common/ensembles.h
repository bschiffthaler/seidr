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

#include <boost/program_options.hpp>
#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif

namespace po = boost::program_options;

template<typename T>
int
check_bootstrap_params(const po::variables_map& vm, T param, seidr_mpi_logger& log)
{
  if (! vm["min-predictor-size"].defaulted()) {
    if (param.predictor_sample_size_min == 0) {
      log << "Looks like you have supplied a fraction to --min-predictor-size"
          << ". This argument takes the minimum count of predictors.\n";
      log.log(LOG_ERR);
      return 1;
    }
  }
  if (! vm["max-predictor-size"].defaulted()) {
    if (param.predictor_sample_size_max == 0) {
      log << "Looks like you have supplied a fraction to --max-predictor-size"
          << ". This argument takes the maximum count of predictors.\n";
      log.log(LOG_ERR);
      return 1;
    }
  }
  if (! vm["min-experiment-size"].defaulted()) {
    if (param.min_sample_size == 0) {
      log << "Looks like you have supplied a fraction to --min-experiment-size"
          << ". This argument takes the minimum count of experiments.\n";
      log.log(LOG_ERR);
      return 1;
    }
  }
  if (! vm["max-experiment-size"].defaulted()) {
    if (param.max_sample_size == 0) {
      log << "Looks like you have supplied a fraction to --max-experiment-size"
          << ". This argument takes the minimum count of experiments.\n";
      log.log(LOG_ERR);
      return 1;
    }
  }
  return 0;
}