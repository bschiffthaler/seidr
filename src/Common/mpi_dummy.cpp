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
#include <BSlogger.hpp>
#include <common.h>
#include <fs.h>
#include <mpi_dummy.h>

#include <sstream>
#include <string>
#include <vector>

// Static members
unsigned int seidr_mpi_logger::_loglevel = LOG_DEFAULT;
bool seidr_mpi_logger::first_instance = true;

seidr_mpi_omp::seidr_mpi_omp(const uint64_t& bs,
                             const arma::mat& data,
                             const std::vector<uint64_t>& indices,
                             const std::vector<std::string>& genes,
                             const std::string& tempdir,
                             const std::string& outfile)
  : _id(0)
  , _procn(0)
  , _pbar(seidr_mpi_progbar<uint64_t>(std::cerr,
                                      indices.size(),
                                      DEFAULT_POLL_INTERVAL,
                                      DEFAULT_BAR_WIDTH,
                                      "genes"))
  , _current_i(0)
  , _bsize(bs)
  , _indices(indices)
  , _data(data)
  , _genes(genes)
  , _outfile(outfile)
  , _tempdir(tempdir)
  , _init_time(0)
  , _queue_fh(nullptr)
{
  // Always take all the work and let OMP handle the scheduling
  _my_indices = _indices;
}

void
seidr_mpi_omp::remove_queue_file()
{}

void
seidr_mpi_omp::get_more_work()
{
  // Since we take all the work at construction time all this needs to do is
  // to set the completion condition
  _my_indices.clear();
}

bool
seidr_mpi_omp::check_logs(
  const std::string& bn) // NOLINT(clang-diagnostic-unused-parameter)
{
  return false;
}

seidr_mpi_logger::seidr_mpi_logger()
  : _rank(0)
{
  _nam = std::string(mpi_get_host());
  if (first_instance) {
    (*this) << "This is seidr v." << _XSTR(SEIDR_VERSION) << '\n';
    this->send(LOG_INFO);
    first_instance = false;
  }
}

seidr_mpi_logger::seidr_mpi_logger(std::string nam)
  : _rank(0)
  , _nam(std::move(nam))
{
  if (first_instance) {
    (*this) << "This is seidr v." << _XSTR(SEIDR_VERSION) << '\n';
    this->send(LOG_INFO);
    first_instance = false;
  }
}

// DEPRECATED: seidr_mpi_logger::send() handles both
// master and slave logging now
void
seidr_mpi_logger::log(unsigned ll)
{
  this->send(ll);
}

void
seidr_mpi_logger::send(unsigned ll)
{
  if (ll > _loglevel) {
    _ss.str(std::string());
    return;
  }

  logger log(std::clog, _nam + ": 0|" + std::to_string(omp_get_thread_num()));
  log(ll) << _ss.str();
  log.flush();
  _ss.str(std::string());
}

void
mpi_sync_tempdir(
  std::string* tempdir) // NOLINT(clang-diagnostic-unused-parameter)
{}
void
mpi_sync_cpr_vector(
  std::vector<uint64_t>* resume) // NOLINT(clang-diagnostic-unused-parameter)
{}

std::string
mpi_get_host()
{
  return std::string("localhost");
}