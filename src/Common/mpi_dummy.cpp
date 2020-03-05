//
// Seidr - Create and operate on gene crowd networks
// Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
//
// This file is part of Seidr.
//
// Seidr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Seidr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
//
#include <BSlogger.h>
#include <common.h>
#include <fs.h>
#include <mpi_dummy.h>

#include <sstream>
#include <string>
#include <vector>

// Static members
unsigned int seidr_mpi_logger::_loglevel = LOG_DEFAULT;

seidr_mpi_omp::seidr_mpi_omp(const uint64_t& bs,
                             const arma::mat& data,
                             const std::vector<uint64_t>& indices,
                             const std::vector<std::string>& genes,
                             const std::string& tempdir,
                             const std::string& outfile) :
  _id(0),
  _procn(0),
  _current_i(0),
  _bsize(bs),
  _indices(indices),
  _data(data),
  _genes(genes),
  _outfile(outfile),
  _tempdir(tempdir),
  _init_time( 0 ),
  _queue_fh(nullptr)
{
  // Always take all the work and let OMP handle the scheduling
  _my_indices = _indices;
}

void seidr_mpi_omp::remove_queue_file() {}

void seidr_mpi_omp::get_more_work() {}

bool seidr_mpi_omp::check_logs(const std::string& bn)
{
  return false;
}

seidr_mpi_logger::seidr_mpi_logger() :
  _rank(0)
{
  _nam = std::string(mpi_get_host());
}

seidr_mpi_logger::seidr_mpi_logger(std::string nam) :
  _rank(0),
  _nam(std::move(nam)) {}

// DEPRECATED: seidr_mpi_logger::send() handles both
// master and slave logging now
void seidr_mpi_logger::log(unsigned ll)
{
  this->send(ll);
}

void seidr_mpi_logger::send(unsigned ll)
{
  if (ll > _loglevel)
  {
    _ss.str(std::string());
    return;
  }

  logger log(std::clog, _nam + ": 0|" +
             std::to_string(omp_get_thread_num()));
  log(ll) << _ss.str();
  log.flush();
  _ss.str(std::string()); 
}

void mpi_sync_tempdir(std::string * tempdir) {}

std::string mpi_get_host()
{
  return std::string("localhost");
}