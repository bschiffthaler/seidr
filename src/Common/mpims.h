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

#include <mpi.h>
#include <armadillo>
#include <vector>
#include <string>
#include <sstream>
#include <fs.h>

#include <BSlogger.h>

//Control messages
#define SEIDR_MPI_READY 0
#define SEIDR_MPI_BUSY 1

#define SEIDR_MPI_ASSIGN 0
#define SEIDR_MPI_EXIT 1
//Tags
#define SEIDR_MPI_RDYS_TAG 1
#define SEIDR_MPI_IDX_TAG 2
#define SEIDR_MPI_CTL_TAG 3
#define SEIDR_MPI_LOG_TAG 4
#define SEIDR_MPI_TEMPDIR_TAG 5

#define SEIDR_MPI_POLL_TIME 100
#define SEIDR_DEF_BATCH 20

class seidr_mpi_logger {
public:
  seidr_mpi_logger();
  template <typename T>
  friend seidr_mpi_logger& operator<<(seidr_mpi_logger& lhs, const T& rhs);
  void send(unsigned ll);
  void log(unsigned ll);
  void set_log_level(unsigned int x) { _loglevel = x; }
  static unsigned int _loglevel;
private:
  std::stringstream _ss;
  int _rank;
};

class seidr_mpi {
public:
  //ctors
  seidr_mpi(const uint64_t& bs, 
            const arma::mat& data, 
            std::vector<unsigned long>& indices,
            const std::vector<std::string>& genes,
            const std::string& tempdir,
            const std::string& outfile);
  //getters
  arma::mat get_data() const {return _data;}
  unsigned long get_current_i() const {return _current_i;}
  //setters
  // void set_data() {_data = m;}
  // void set_indices() {_indices = u;}
  // void set_genes() {_genes = g;}
  // void set_outfilebase() {_outfilebase = b;}
  // void set_tempdir() {_tempdir = b;}
  // void set_outfile() {_outfile = b;}
  //main function
  void exec();
  virtual void entrypoint();
  virtual void finalize();
protected:
  //MPI constants
  int _id;
  int _procn;
  //MPI ready states
  int _readystate;
  std::vector<int> _idle;
  void poll_idle_proc();
  bool nb_poll_idle_proc(); //non-blocking
  void announce_ready();
  int next_idle();
  int nb_next_idle();
  //MPI logging
  void send_log(std::string s);
  bool check_logs();
  //
  unsigned long _current_i;
  unsigned long _bsize;
  std::vector<unsigned long>& _indices;
  const arma::mat& _data;
  const std::vector<std::string>& _genes;
  const std::string& _outfile;
  const std::string& _tempdir;
  //misc
  double _init_time;
};

template <typename T>
seidr_mpi_logger&
operator<<(seidr_mpi_logger& lhs, const T& rhs)
{
  lhs._ss << rhs;
  return lhs;
}

void mpi_sync_tempdir(std::string * tempdir);
