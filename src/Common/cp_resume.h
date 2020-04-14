#pragma once

#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif
#include <common.h>
#include <fs.h>

#include <algorithm>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = boost::filesystem;

enum cpr_format_t
{
  CPR_M,
  CPR_LM
};

template<typename T>
class cp_resume
{
public:
  cp_resume(T& param, cpr_format_t format)
  {
    _param = param;
    _format = format;
  }
  void print();
  void save(const std::string& file, const T& param);
  void load(T& param, const std::string& file);
  void resume();
  std::vector<uint64_t> get_good_idx() { return _good_idx; }

private:
  bool _valid(const uint64_t& nf, const uint64_t& ng, const uint64_t& i);
  T _param;
  cpr_format_t _format;
  std::vector<uint64_t> _good_idx;
  std::string _resume_file;
};

template<typename T>
bool
cp_resume<T>::_valid(const uint64_t& nf, const uint64_t& ng, const uint64_t& i)
{
  if (_format == CPR_M) {
    return nf == ng;
  } else if (_format == CPR_LM) {
    return nf == i;
  } else {
    _BUG("Control reaches unexpected code block.");
    return false;
  }
}

template<typename T>
void
cp_resume<T>::print()
{
  std::stringstream ss;
  boost::archive::xml_oarchive oa(ss);
  oa << BOOST_SERIALIZATION_NVP(_param);
  std::cerr << ss.str();
}

template<typename T>
void
cp_resume<T>::save(const std::string& file, const T& param)
{
  std::ofstream ofs(file.c_str());

  if (!ofs.good()) {
    throw std::runtime_error("Could not open " + file + " for writing.");
  }

  boost::archive::xml_oarchive oa(ofs);
  oa << BOOST_SERIALIZATION_NVP(param);

  ofs.flush();
}

template<typename T>
void
cp_resume<T>::load(T& param, const std::string& file)
{
  seidr_mpi_logger log;
  _resume_file = file;

  std::ifstream ifs(file.c_str());

  if (!ifs.good()) {
    throw std::runtime_error("Could not open " + file + " for reading.");
  }

  boost::archive::xml_iarchive ia(ifs);
  ia >> BOOST_SERIALIZATION_NVP(param);

  param.resuming = true;
  _param = param;
}

template<typename T>
void
cp_resume<T>::resume()
{
  seidr_mpi_logger log("CPR@" + mpi_get_host());
  log << "Checking completed files...\n";
  log.log(LOG_INFO);
  std::vector<fs::path> files;
  fs::path p_tmp(_param.tempdir);

  // Collect all files in temp directory that lloosely follow naming convention
  for (auto it = fs::directory_iterator(p_tmp); it != fs::directory_iterator();
       it++) {
    std::string pstring = it->path().string();
    if (pstring.find(".") != std::string::npos) {
      log << "Ignoring unexpected file: " << pstring << '\n';
      log.log(LOG_WARN);
    } else if (fs::is_regular_file(it->path())) {
      files.push_back((*it).path());
    }
  }

  std::vector<std::string> genes = read_genes(_param.gene_file);
  std::string tmpf = tempfile(_param.tempdir);
  std::ofstream tmpf_ofs(tmpf.c_str());

  if (!tmpf_ofs.good()) {
    throw std::runtime_error("Unable to create CPR tempfile: " + tmpf);
  }

  for (const fs::path& p : files) {
    uint64_t i = 0;
    uint64_t cur_idx = 0;
    uint64_t cur_nf = 0;
    std::string line;
    std::ifstream ifs(p.string().c_str());
    while (std::getline(ifs, line)) {
      if (i % 2 == 0) {
        cur_idx = std::stoul(line);
      } else {
        std::string fields;
        std::stringstream ss(line);
        std::vector<double> values;
        while (std::getline(ss, fields, '\t')) {
          try {
            values.push_back(std::stod(fields));
          } catch (std::exception& e) {
            continue;
          }
          cur_nf++;
        }
        if (_valid(cur_nf, genes.size(), cur_idx)) {
          _good_idx.push_back(cur_idx);
          tmpf_ofs << cur_idx << '\n';
          for (uint64_t i2 = 0; i2 < values.size(); i2++) {
            tmpf_ofs << values[i2];
            tmpf_ofs << ((i2 == (values.size() - 1)) ? '\n' : '\t');
          }
          cur_idx = 0;
          cur_nf = 0;
        }
      }
      i++;
    }
    fs::remove(p);
  }
  tmpf_ofs.close();
  std::sort(_good_idx.begin(), _good_idx.end());
  log << "Rescued " << _good_idx.size() << " genes\n";
  log.log(LOG_INFO);
}