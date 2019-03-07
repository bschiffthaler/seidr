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

#include <mpi.h>
#include <armadillo>
#include <vector>
#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>

#include <BSlogger.h>
#include <mpims.h>

using boost::lexical_cast;
unsigned int seidr_mpi_logger::_loglevel = LOG_DEFAULT;

seidr_mpi_logger::seidr_mpi_logger() :
  _rank(0)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &_rank); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
}

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

  if (_rank == 0)
  {
    logger log(std::clog, "Master");
    log(ll) << _ss.str();
    log.flush();
    _ss.str(std::string());
  }
  else
  {
    std::string l;
    switch (ll)
    {
    case LOG_DEBUG:
      l = "D"; break;
    case LOG_WARN:
      l = "W"; break;
    case LOG_ERR:
      l = "E"; break;
    case LOG_INFO:
      l = "I"; break;
    case LOG_TIME:
      l = "T"; break;
    default:
      l = "x"; break;
    }
    std::string s = l + _ss.str();
    MPI_Send(&s[0], s.size(), MPI_CHAR, 0, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
             SEIDR_MPI_LOG_TAG, MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    _ss.str(std::string());
  }
}

// seidr_mpi::seidr_mpi() :
//   _id(0),
//   _procn(0),
//   _readystate(SEIDR_MPI_READY),
//   _current_i(0),
//   _bsize(SEIDR_DEF_BATCH),
//   _init_time( MPI_Wtime() )
// {
//   MPI_Comm_rank(MPI_COMM_WORLD, &_id); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
//   MPI_Comm_size(MPI_COMM_WORLD, &_procn); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
//   if (_id == 0) //master
//   {}
//   else //slave
//   {
//     announce_ready();
//   }
// }

seidr_mpi::seidr_mpi(const uint64_t& bs,
                     const arma::mat& data,
                     std::vector<unsigned long>& indices,
                     const std::vector<std::string>& genes,
                     const std::string& tempdir,
                     const std::string& outfile) :
  _id(0),
  _procn(0),
  _readystate(SEIDR_MPI_READY),
  _current_i(0),
  _bsize(bs),
  _indices(indices),
  _data(data),
  _genes(genes),
  _outfile(outfile),
  _tempdir(tempdir),
  _init_time( MPI_Wtime() )
{
  MPI_Comm_rank(MPI_COMM_WORLD, &_id); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  MPI_Comm_size(MPI_COMM_WORLD, &_procn); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  if (_id == 0) //master
  {}
  else //slave
  {
    announce_ready();
  }
}

void seidr_mpi::announce_ready()
{
  MPI_Send(&_readystate, 1, MPI_INT, 0, SEIDR_MPI_RDYS_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
           MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
}

void seidr_mpi::poll_idle_proc()
{
  if (_procn != 1 && _id == 0)
  {
    MPI_Status status;
    int rdystate;
    MPI_Recv(&rdystate, 1, MPI_INT, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
             MPI_ANY_SOURCE, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
             SEIDR_MPI_RDYS_TAG,
             MPI_COMM_WORLD, &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    _idle.push_back(status.MPI_SOURCE);
  }
}

bool seidr_mpi::nb_poll_idle_proc()
{
  bool ret = false;
  if (_procn == 1)
  {
    ret = false;
  }
  else if (_id == 0)
  {
    MPI_Status probe_status;
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, SEIDR_MPI_RDYS_TAG,
               MPI_COMM_WORLD, &flag, &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    if (flag == 1)
    {
      MPI_Status status;
      int rdystate;
      MPI_Recv(&rdystate, 1, MPI_INT, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               MPI_ANY_SOURCE, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               SEIDR_MPI_RDYS_TAG,
               MPI_COMM_WORLD, &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      if (rdystate == SEIDR_MPI_READY)
      {
        _idle.push_back(status.MPI_SOURCE);
      }
      ret = true;
    }
  }
  return ret;
}

int seidr_mpi::next_idle()
{
  int ret = 0;
  if (_procn == 1)
  {
    ret = 0;
  }
  else if (_idle.empty())
  {
    poll_idle_proc();
    int it = _idle[_idle.size() - 1];
    _idle.pop_back();
    ret = it;
  }
  else
  {
    int it = _idle[_idle.size() - 1];
    _idle.pop_back();
    ret = it;
  }
  return ret;
}

int seidr_mpi::nb_next_idle()
{
  if (_procn == 1)
  {
    return 0;
  }

  while (_idle.empty())
  {
    nb_poll_idle_proc();
  }
  int it = _idle[_idle.size() - 1];
  _idle.pop_back();
  return it;
}


bool seidr_mpi::check_logs()
{
  if (_id == 0)
  {
    MPI_Status probe_status;
    int flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, SEIDR_MPI_LOG_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               MPI_COMM_WORLD, &flag, &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    if (flag == 1)
    {
      std::string s;
      int size;
      MPI_Status status;
      int source = probe_status.MPI_SOURCE;
      MPI_Get_count(&probe_status, MPI_CHAR, &size); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      s.resize(size);
      MPI_Recv(&s[0], size, MPI_CHAR, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               source,
               SEIDR_MPI_LOG_TAG,
               MPI_COMM_WORLD, &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      unsigned l;
      switch (s[0])
      {
      case 'D':
        l = LOG_DEBUG; break;
      case 'W':
        l = LOG_WARN; break;
      case 'E':
        l = LOG_ERR; break;
      case 'I':
        l = LOG_INFO; break;
      case 'T':
        l = LOG_TIME; break;
      default:
        l = LOG_DEFAULT; break;
      }
      s.erase(0, 1);
      std::string loc =
        (source == 0 ?
         "Master" :
         ("Thread " + std::to_string(source) ) );
      logger log(std::clog, loc);
      log(l) << s;

      flag = 0;
      MPI_Iprobe(MPI_ANY_SOURCE, SEIDR_MPI_LOG_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
                 MPI_COMM_WORLD, &flag, &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      return true;
    }
  }
  return false;
}
void seidr_mpi::exec()
{
  if (_id == 0 && _procn > 1) //master
  {
    seidr_mpi_logger log;
    double start_time = MPI_Wtime();
    while (_current_i < _indices.size())
    {
      // Process logs
      check_logs();

      int curr_thread;
      if (! _idle.empty())
      {
        curr_thread = next_idle();
      }
      else
      {
        poll_idle_proc();
        while (_idle.empty())
        {
          poll_idle_proc();
          check_logs();
        }
        curr_thread = next_idle();
      }

      int ctl = SEIDR_MPI_ASSIGN;
      MPI_Send(&ctl,
               1,
               MPI_INT, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               curr_thread,
               SEIDR_MPI_CTL_TAG,
               MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      if ((_indices.size() - _current_i) < _bsize)
      {
        int count = (_indices.size() - _current_i);
        log << "Indices " << _current_i << "-"
            << (_current_i + count) - 1 << " assigned to thread "
            << curr_thread << '\n';
        log.send(LOG_INFO);
        MPI_Send(&_indices[_current_i],
                 count,
                 MPI_UNSIGNED_LONG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
                 curr_thread,
                 SEIDR_MPI_IDX_TAG,
                 MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
        _current_i += count;
      }
      else
      {
        log << "Indices " << _current_i << "-"
            << (_current_i + _bsize) - 1
            << " assigned to thread "
            << curr_thread << '\n';
        log.send(LOG_INFO);
        MPI_Send(&_indices[_current_i],
                 _bsize,
                 MPI_UNSIGNED_LONG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
                 curr_thread,
                 SEIDR_MPI_IDX_TAG,
                 MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
        _current_i += _bsize;
      }
    }
    // Finish up
    while (lexical_cast<int>(_idle.size()) < _procn - 1)
    {
      nb_poll_idle_proc();
      check_logs();
    }
    // Set exit flags on idle threads
    for (int i = 0; i < _procn - 1; i++)
    {
      int curr_thread = next_idle();
      int ctl = SEIDR_MPI_EXIT;
      MPI_Send(&ctl,
               1,
               MPI_INT, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               curr_thread,
               SEIDR_MPI_CTL_TAG,
               MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    }
    double finish_time = MPI_Wtime();
    log << "Finished computation in " << finish_time - start_time
        << " seconds\n";
    log.send(LOG_INFO);
    finalize();
    MPI_Finalize();
    finish_time = MPI_Wtime();
    log << "Finished cleanup in " << finish_time - start_time
        << " seconds\n";
    log.send(LOG_INFO);
  }
  else if (_procn > 1) // slave
  {
    bool exit = false;
    while (!exit)
    {
      _indices.clear();
      // Get control tag (more data or shutdown)
      int ctl;
      MPI_Status status;
      MPI_Recv(&ctl, 1, MPI_INT, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
               0,
               SEIDR_MPI_CTL_TAG,
               MPI_COMM_WORLD, &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      // Assign work to thread
      if (ctl == SEIDR_MPI_ASSIGN)
      {
        MPI_Status probe;
        MPI_Probe(0, SEIDR_MPI_IDX_TAG, MPI_COMM_WORLD, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
                  &probe);
        int bufsize;
        MPI_Get_count(&probe, MPI_UNSIGNED_LONG, &bufsize); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
        _indices.resize(bufsize);
        MPI_Recv(&_indices[0], bufsize, MPI_UNSIGNED_LONG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
                 0,
                 SEIDR_MPI_IDX_TAG,
                 MPI_COMM_WORLD, &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
        entrypoint();
      }
      // Exit loop on exit message
      else if (ctl == SEIDR_MPI_EXIT)
      {
        exit = true;
      }
    }
    MPI_Finalize();
  }
  else // We are running single threaded
  {
    seidr_mpi_logger log;
    double start_time = MPI_Wtime();
    entrypoint();
    finalize();
    double finish_time = MPI_Wtime();
    log << "Finished. Took " << finish_time - start_time
        << " seconds\n";
    log.send(LOG_INFO);
    MPI_Finalize();
  }
}

void seidr_mpi::entrypoint()
{
  announce_ready();
}

void seidr_mpi::finalize()
{}

void mpi_sync_tempdir(std::string * tempdir)
{
  int ntasks;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)

  if (rank == 0)
  {
    if (ntasks > 1)
    {
#ifdef DEBUG
      log << "Informing threads of tempdir: " << (*tempdir) << '\n';
      log.log(LOG_DEBUG);
#endif
      for (int i = 1; i < ntasks; i++)
      {
        MPI_Send(&(*tempdir)[0], tempdir->size(), MPI_CHAR, i, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
                 SEIDR_MPI_TEMPDIR_TAG, MPI_COMM_WORLD); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
      }
    }
  }
  else
  {
    MPI_Status probe_status;
    MPI_Probe(MPI_ANY_SOURCE, SEIDR_MPI_TEMPDIR_TAG, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
              MPI_COMM_WORLD, &probe_status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    int size;
    MPI_Status status;
    int source = probe_status.MPI_SOURCE;
    MPI_Get_count(&probe_status, MPI_CHAR, &size); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    tempdir->resize(size);
    MPI_Recv(&(*tempdir)[0], size, MPI_CHAR, // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
             source,
             SEIDR_MPI_TEMPDIR_TAG,
             MPI_COMM_WORLD, &status); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
#ifdef DEBUG
    log << "Got informed of tempdir: " << (*tempdir) << '\n';
    log.log(LOG_DEBUG);
#endif
  }
}
