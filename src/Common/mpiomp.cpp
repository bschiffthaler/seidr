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
#include <mpi.h>
#include <mpiomp.h>

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
  , _my_indices(std::vector<uint64_t>())
  , _data(data)
  , _genes(genes)
  , _outfile(outfile)
  , _tempdir(tempdir)
  , _init_time(MPI_Wtime())
  , _queue_file(std::string())
  , _queue_fh(nullptr)
{
  MPI_Comm_rank(MPI_COMM_WORLD,
                &_id); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  MPI_Comm_size(MPI_COMM_WORLD,
                &_procn); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  // Make sure everyoneuses the same queue file
  std::string queue_file = tempfile(tempdir) + ".q";
  mpi_sync_tempdir(&queue_file);
  _queue_file = queue_file;
  MPI_Barrier(MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  if (_id == 0) {
    MPI_File ftmp;
    MPI_File_open(
      MPI_COMM_SELF, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      _queue_file.c_str(),
      MPI_MODE_CREATE | MPI_MODE_WRONLY, // NOLINT
      MPI_INFO_NULL,
      &ftmp); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    MPI_Status status;
    MPI_File_write(ftmp,
                   &_indices[0],
                   _indices.size(),
                   MPI_UINT64_T,
                   &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    MPI_File_close(&ftmp);
  }
  // Everyone needs to wait until queue is set up
  MPI_Barrier(MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  MPI_File_open(
    MPI_COMM_WORLD, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    _queue_file.c_str(),
    MPI_MODE_RDONLY, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    MPI_INFO_NULL,
    &_queue_fh); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  get_more_work();
}

void
seidr_mpi_omp::remove_queue_file()
{
  // Absolutely all MPI tasks need to finish before this point
  MPI_Barrier(MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  MPI_File_close(&_queue_fh);
  MPI_Barrier(MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  if (_id == 0) {
    remove(_queue_file.c_str());
  }
}

void
seidr_mpi_omp::get_more_work()
{
  _my_indices.resize(_bsize);
  MPI_Status status;

  MPI_File_read_shared(
    _queue_fh,
    &_my_indices[0],
    _bsize,
    MPI_UINT64_T,
    &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  int work;
  MPI_Get_count(&status,
                MPI_UINT64_T,
                &work); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  _my_indices.resize(work);
}

bool
seidr_mpi_omp::check_logs(const std::string& bn)
{
  if (_id == 0) {
    MPI_Status probe_status;
    int flag = 0;
    MPI_Iprobe(
      MPI_ANY_SOURCE,
      SEIDR_MPI_LOG_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      MPI_COMM_WORLD,
      &flag,
      &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    if (flag == 1) {
      std::string s;
      int size;
      MPI_Status status;
      int source = probe_status.MPI_SOURCE;
      MPI_Get_count(&probe_status,
                    MPI_CHAR,
                    &size); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      s.resize(size);
      MPI_Recv(&s[0],
               size,
               MPI_CHAR, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               source,
               SEIDR_MPI_LOG_TAG,
               MPI_COMM_WORLD,
               &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      unsigned l;
      switch (s[0]) {
        case 'D':
          l = LOG_DEBUG;
          break;
        case 'W':
          l = LOG_WARN;
          break;
        case 'E':
          l = LOG_ERR;
          break;
        case 'I':
          l = LOG_INFO;
          break;
        case 'T':
          l = LOG_TIME;
          break;
        default:
          l = LOG_DEFAULT;
          break;
      }
      s.erase(0, 1);
      std::string loc = bn + ": " + std::to_string(source) + "|" +
                        std::to_string(omp_get_thread_num());
      logger log(std::clog, loc);
      log(l) << s;

      flag = 0;
      MPI_Iprobe(
        MPI_ANY_SOURCE,
        SEIDR_MPI_LOG_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
        MPI_COMM_WORLD,
        &flag,
        &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      return true;
    }
  }
  return false;
}

seidr_mpi_logger::seidr_mpi_logger()
  : _rank(0)
{
  _nam = std::string(mpi_get_host());
  MPI_Comm_rank(MPI_COMM_WORLD,
                &_rank); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  (*this) << "This is seidr v." << _XSTR(SEIDR_VERSION) << '\n';
  this->send(LOG_INFO); 
}

seidr_mpi_logger::seidr_mpi_logger(const std::string & nam)
  : _rank(0)
  , _nam(std::move(nam))
{
  MPI_Comm_rank(MPI_COMM_WORLD,
                &_rank); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  (*this) << "This is seidr v." << _XSTR(SEIDR_VERSION) << '\n';
  this->send(LOG_INFO); 
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

  if (_rank == 0) {
    logger log(std::clog, _nam + ": 0|" + std::to_string(omp_get_thread_num()));
    log(ll) << _ss.str();
    log.flush();
    _ss.str(std::string());
  } else {
    std::string l;
    switch (ll) {
      case LOG_DEBUG:
        l = "D";
        break;
      case LOG_WARN:
        l = "W";
        break;
      case LOG_ERR:
        l = "E";
        break;
      case LOG_INFO:
        l = "I";
        break;
      case LOG_TIME:
        l = "T";
        break;
      default:
        l = "x";
        break;
    }
    std::string s = l + _ss.str();
    MPI_Send(&s[0],
             s.size(),
             MPI_CHAR,
             0, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
             SEIDR_MPI_LOG_TAG,
             MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    _ss.str(std::string());
  }
}

void
mpi_sync_tempdir(std::string* tempdir)
{
  int ntasks;
  MPI_Comm_size(MPI_COMM_WORLD,
                &ntasks); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,
                &rank); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

#ifdef DEBUG
  seidr_mpi_logger log("common@" + mpi_get_host());
#endif

  if (rank == 0) {
    if (ntasks > 1) {
#ifdef DEBUG
      log << "Informing threads of tempdir: " << (*tempdir) << '\n';
      log.send(LOG_DEBUG);
#endif
      for (int i = 1; i < ntasks; i++) {
        MPI_Send(
          &(*tempdir)[0],
          tempdir->size(),
          MPI_CHAR,
          i, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
          SEIDR_MPI_TEMPDIR_TAG,
          MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      }
    }
  } else {
    MPI_Status probe_status;
    MPI_Probe(
      MPI_ANY_SOURCE,
      SEIDR_MPI_TEMPDIR_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      MPI_COMM_WORLD,
      &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    int size;
    MPI_Status status;
    int source = probe_status.MPI_SOURCE;
    MPI_Get_count(&probe_status,
                  MPI_CHAR,
                  &size); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    tempdir->resize(size);
    MPI_Recv(&(*tempdir)[0],
             size,
             MPI_CHAR, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
             source,
             SEIDR_MPI_TEMPDIR_TAG,
             MPI_COMM_WORLD,
             &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
#ifdef DEBUG
    log << "Got informed of tempdir: " << (*tempdir) << '\n';
    log.log(LOG_DEBUG);
#endif
  }
}

void
mpi_sync_cpr_vector(std::vector<uint64_t>* resume)
{
  int ntasks;
  MPI_Comm_size(MPI_COMM_WORLD,
                &ntasks); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,
                &rank); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

#ifdef DEBUG
  seidr_mpi_logger log("common@" + mpi_get_host());
#endif

  if (rank == 0) {
    if (ntasks > 1) {
      for (int i = 1; i < ntasks; i++) {
        MPI_Send(
          &(*resume)[0],
          resume->size(),
          MPI_UNSIGNED_LONG,
          i, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
          SEIDR_MPI_CPR_TAG,
          MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      }
    }
  } else {
    MPI_Status probe_status;
    MPI_Probe(
      MPI_ANY_SOURCE,
      SEIDR_MPI_CPR_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      MPI_COMM_WORLD,
      &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    int size;
    MPI_Status status;
    int source = probe_status.MPI_SOURCE;
    MPI_Get_count(&probe_status,
                  MPI_UNSIGNED_LONG,
                  &size); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    resume->resize(size);
    MPI_Recv(
      &(*resume)[0],
      size,
      MPI_UNSIGNED_LONG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      source,
      SEIDR_MPI_CPR_TAG,
      MPI_COMM_WORLD,
      &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
#ifdef DEBUG
    log << "Got informed of tempdir: " << (*tempdir) << '\n';
    log.log(LOG_DEBUG);
#endif
  }
}

std::string
mpi_get_host()
{
  std::string host;
  host.resize(
    MPI_MAX_PROCESSOR_NAME); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  int hostlen;
  MPI_Get_processor_name(&host[0], &hostlen);
  host.resize(hostlen);
  return host;
}