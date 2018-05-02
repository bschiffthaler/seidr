#pragma once

#include <mpi.h>
#include <armadillo>
#include <vector>
#include <string>
#include <sstream>

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

#define SEIDR_MPI_POLL_TIME 100
#define SEIDR_DEF_BATCH 20

class seidr_mpi_logger {
 public:
  seidr_mpi_logger();
  template <typename T>
    friend seidr_mpi_logger& operator<<(seidr_mpi_logger& lhs, const T& rhs);
  void send(unsigned ll);
  void log(unsigned ll);
 private:
  std::stringstream _ss;
  int _rank;
  unsigned int _loglevel;
};

class seidr_mpi {
public:
  //ctors
  seidr_mpi();
  seidr_mpi(unsigned long bs);
  //getters
  arma::mat get_data() const {return _data;}
  unsigned long get_current_i() const {return _current_i;}
  //setters
  void set_data(arma::mat m){_data = m;}
  void set_indices(std::vector<unsigned long> u){_indices = u;}
  void set_genes(std::vector<std::string> g){_genes = g;}
  void set_outfilebase(std::string b) {_outfilebase = b;}
  void set_outfile(std::string b) {_outfile = b;}
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
  std::vector<unsigned long> _indices;
  arma::mat _data;
  std::vector<std::string> _genes;
  std::string _outfile;
  std::string _outfilebase;
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
