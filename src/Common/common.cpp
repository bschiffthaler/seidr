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

// Seidr
#include <Serialize.h>
#include <common.h>
#include <fs.h>
#ifdef SEIDR_WITH_MPI
#include <mpiomp.h>
#else
#include <mpi_dummy.h>
#endif
// External
#include <cerrno>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using file_index_t = std::pair<std::streampos, std::string>;
using result_file_map_t = std::map<seidr_uword_t, file_index_t>;

arma::uvec
get_i(arma::uword ind, size_t s)
{
  arma::uvec res(s - 1);
  arma::uword ri = 0;
  for (arma::uword f = 0; f < s; f++) {
    if (f < ind) {
      res(ri++) = f;
    }
    if (f > ind) {
      res(ri++) = f;
    }
  }
  return res;
}

std::vector<std::string>
read_genes(const std::string& input, char row_delim, char field_delim)
{
  std::vector<std::string> res;
  const char rd = row_delim;
  const char fd = field_delim;

  std::ifstream inf;
  inf.open(input.c_str());
  for (std::string row; std::getline(inf, row, rd);) {
    std::istringstream ss(row);

    for (std::string field; std::getline(ss, field, fd);) {
      res.push_back(field);
    }
  }

  inf.close();

  return res;
}

bool
is_gzip(const std::string& input)
{
  std::ifstream _input(input.c_str(), std::ios::in | std::ios::binary);
  bool gzip = false;
  char byte1 = _input.get();
  char byte2 = _input.get();
  if ((byte1 == '\x1F') && (byte2 == '\x8B')) {
    gzip = true;
  }
  _input.close();
  return gzip;
}

// Almost equal a la Bruce Dawson
bool
almost_equal(seidr_score_t A, seidr_score_t B)
{
  // Calculate the difference.
  seidr_score_t diff = fabs(A - B);
  A = fabs(A);
  B = fabs(B);
  // Find the largest
  seidr_score_t largest = (B > A) ? B : A;

  return (diff <= largest * SFLT_EPSILON);
}

seidr_score_t
unity_stand(seidr_score_t xmin, seidr_score_t xmax, seidr_score_t xi)
{
  seidr_score_t x = 1 - ((xi - xmin) / (xmax - xmin));
  return x;
}

void
scale(arma::mat& x)
{
  x.each_col([](arma::vec& v) {
    double m = arma::mean(v);
    double sd = arma::stddev(v);
    v -= m;
    v /= sd;
  });
}

std::vector<std::string>
tokenize_delim(const std::string& nodes, const std::string& delim)
{
  std::vector<std::string> nodelist;
  boost::char_separator<char> sep(delim.c_str());
  tokenizer tokens(nodes, sep);
  for (auto tok = tokens.begin(); tok != tokens.end(); tok++) {
    nodelist.push_back(*tok);
  }
  return nodelist;
}

void
merge_files(const std::string& outfile,
            std::string tempdir,
            bool targeted,
            int id,
            const std::vector<std::string>& genes)
{
  if (id == 0) {
    result_file_map_t rmap;
    seidr_mpi_logger log("merge_files@" + mpi_get_host());
    log << "Looking for files to merge...\n";
    log.send(LOG_INFO);
    std::vector<fs::path> files;
    fs::path p_tmp(tempdir);
    // Collect all files in temp directory that lloosely follow naming
    // convention
    for (auto it = fs::directory_iterator(p_tmp);
         it != fs::directory_iterator();
         it++) {
      std::string pstring = it->path().string();
      if (pstring.find(".") != std::string::npos) {
        log << "Ignoring unexpected file: " << pstring << '\n';
        log.log(LOG_WARN);
      } else if (fs::is_regular_file(it->path())) {
        files.push_back((*it).path());
      }
    }
    std::ofstream ofs(outfile);

    // Collect file offset and gene ID of all targets in temp file
    for (fs::path& p : files) {
      std::ifstream ifs(p.string().c_str());
      std::string l;
      while (std::getline(ifs, l)) {
        seidr_uword_t gene_index = std::stoul(l);
        std::streampos g = ifs.tellg();
        rmap[gene_index] = file_index_t(g, p.string());
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      ifs.close();
    }

    log << "Merging " << rmap.size() << " genes from " << files.size()
        << " temporary files\n";
    log.send(LOG_INFO);

    std::ifstream ifs;
    std::string file_path;
    for (auto it = rmap.begin(); it != rmap.end(); it++) {
      if (file_path.empty()) // First iteration
      {
        file_path = it->second.second;
        ifs = std::ifstream(it->second.second.c_str());
      } else if (file_path != it->second.second) // New file
      {
        ifs.close();
        file_path = it->second.second;
        ifs = std::ifstream(it->second.second.c_str());
      }
      std::streampos gene_index = it->second.first;

      ifs.seekg(gene_index);
      std::string l;
      std::getline(ifs, l);

      if (targeted) {
        std::stringstream ss(l);
        std::string token;
        seidr_uword_t j = 0;
        seidr_uword_t i = it->first;
        while (ss >> token) {
          if (i != j) {
            ofs << genes[i] << '\t' << genes[j] << '\t' << token << '\n';
          }
          j++;
        }
      } else {
        ofs << l << '\n';
      }
    }
    ofs.close();
    ifs.close();
    for (auto it = fs::directory_iterator(p_tmp);
         it != fs::directory_iterator();
         it++) {
      remove(it->path().string());
    }
    remove(tempdir, true);
  }
}

bool
any_const_expr(arma::mat& inp)
{
  arma::mat v = arma::var(inp);
  for (arma::uword i = 0; i < v.n_elem; i++) {
    if (almost_equal(v(0, i), 0)) {
      return true;
    }
  }
  return false;
}

void
verify_matrix(arma::mat& inp, std::vector<std::string>& genes)
{
  if (genes.size() != inp.n_cols) {
    throw std::runtime_error("There must be as many gene names as columns "
                             "in the expression matrix.");
  }
  verify_matrix(inp);
}

void
verify_matrix(arma::mat& inp)
{
  if (!inp.is_finite()) {
    throw std::runtime_error("Not all elements in input matrix are finite.");
  }
  if (any_const_expr(inp)) {
    throw std::runtime_error("Constant values detected in at least one column"
                             ". Please filter your input to contain only "
                             "columns with non-zero variance.");
  }
}

void
assert_mutually_exclusive(const po::variables_map& vm,
                          const std::vector<std::string> targets)
{
  uint64_t count = 0;
  for (const auto& t : targets) {
    if (vm.count(t) > 0) {
      if (!vm[t].defaulted()) {
        count++;
      }
    }
  }
  if (count > 1) {
    throw std::runtime_error("Arguments " + str_join(targets, ",") +
                             " are mutually exclusive");
  }
}

std::string
str_join(const std::vector<std::string>& source, const std::string& delim)
{
  std::string ret;
  for (uint64_t i = 0; i < source.size(); i++) {
    ret += source[i];
    if (i < (source.size() - 1)) {
      ret += delim;
    }
  }
  return ret;
}

void
make_tpos(uint32_t& tpos, const SeidrFileHeader& h)
{
  if (tpos == 0) {
    tpos = h.attr.nalgs - 1;
  } else {
    tpos--;
  }
  assert_in_range<uint32_t>(tpos, 0, h.attr.nalgs - 1);
}

uint64_t
get_mpi_nthread()
{
  int procn = 0;
#ifdef SEIDR_WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &procn);
#endif
  if (procn < 0) {
    throw std::runtime_error("Number of MPI threads could not be"
                             " reliably determined");
  }
  uint64_t ret = static_cast<uint64_t>(procn);
  return ret;
}

uint64_t
guess_batch_size(uint64_t const& set_size, uint64_t const& task_n)
{
  if (task_n == 0) {
    return set_size;
  }
  if (set_size % task_n == 0) {
    return set_size / task_n;
  } else {
    return (set_size / task_n) + 1;
  }
}