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

#include <BSlogger.hpp>
#include <common.h>
#include <fs.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

#define SEIDR_MPI_LOG_TAG 0
#define SEIDR_MPI_TEMPDIR_TAG 1
#define SEIDR_MPI_CPR_TAG 2
#define SEIDR_MPI_PROGBAR_TAG 3

constexpr uint64_t DEFAULT_POLL_INTERVAL = 1000;
constexpr uint64_t DEFAULT_BAR_WIDTH = 25;

template<typename T>
class seidr_mpi_progbar
{
public:
  seidr_mpi_progbar(std::ostream& f,
                    T max,
                    uint64_t poll_interval = DEFAULT_POLL_INTERVAL,
                    uint64_t width = DEFAULT_BAR_WIDTH,
                    std::string unit = "units")
    : _pbar(progbar_fancy<T>(f, max, poll_interval, width, unit))
    , _rank(0)
  {}
  seidr_mpi_progbar& operator++()
  {
    _pbar++;
    return *this;
  }
  void finalize() { _pbar.finalize(); }
  seidr_mpi_progbar operator++(int)
  {
    const seidr_mpi_progbar<T> copy(*this);
    _pbar++;
    return copy;
  }

private:
  progbar_fancy<T> _pbar;
  int _rank;
  std::string _nam;
  std::string _host;
};

class seidr_mpi_omp
{
public:
  // ctors
  seidr_mpi_omp(const uint64_t& bs,
                const arma::mat& data,
                const std::vector<uint64_t>& indices,
                const std::vector<std::string>& genes,
                const std::string& tempdir,
                const std::string& outfile);
  // getters
  arma::mat get_data() const { return _data; }
  uint64_t get_current_i() const { return _current_i; }
  int rank() { return 0; }
  bool check_logs(const std::string& bn);
  void get_more_work();
  void remove_queue_file();
  void finalize_pbar() { _pbar.finalize(); }
  void increment_pbar() { _pbar++; }

protected:
  // MPI constants
  int _id;
  int _procn;
  seidr_mpi_progbar<uint64_t> _pbar;
  //
  uint64_t _current_i;
  uint64_t _bsize;
  const std::vector<uint64_t>& _indices;
  std::vector<uint64_t> _my_indices;
  const arma::mat& _data;
  const std::vector<std::string>& _genes;
  const std::string& _outfile;
  const std::string& _tempdir;
  // misc
  double _init_time;
  std::string _queue_file; // can't be const as it needs to sync
  int* _queue_fh;
};

class seidr_mpi_logger
{
public:
  seidr_mpi_logger();
  explicit seidr_mpi_logger(std::string nam);
  template<typename T>
  friend seidr_mpi_logger& operator<<(seidr_mpi_logger& lhs, const T& rhs);
  void send(unsigned ll);
  void log(unsigned ll);
  void static set_log_level(unsigned int x) { _loglevel = x; }
  static unsigned int _loglevel;

private:
  std::stringstream _ss;
  int _rank;
  std::string _nam;
  std::string _host;
};

template<typename T>
seidr_mpi_logger&
operator<<(seidr_mpi_logger& lhs, const T& rhs)
{
  lhs._ss << rhs;
  return lhs;
}

void
mpi_sync_tempdir(std::string* tempdir);
void
mpi_sync_cpr_vector(std::vector<uint64_t>* resume);

std::string
mpi_get_host();