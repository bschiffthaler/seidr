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

#include <armadillo>
#include <bgzf.h>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#define HEADER_DIM_SIZE 128
#define HEADER_NALGS_SIZE 16
#define HEADER_VERSION_SIZE 128
#define HEADER_CMD_SIZE 65536
#define HEADER_NODE_SIZE 256
#define HEADER_ALG_SIZE 256
#define HEADER_SUPP_SIZE 128
#define EDGE_DESCRIPTOR_SIZE 64
#define EDGE_DIRECTION_SIZE 8
#define EDGE_SCORE_SIZE 128
#define EDGE_PRESENCE_SIZE 8
#define EDGE_SUPP_SIZE 8

// 00000000
// -------1 Edge exists
// ------1- Edge is directed A->B
// -----1-- Edge is directed A<-B

inline bool
EDGE_EXISTS(const uint8_t& E)
{
  return (E & 0x1) != 0;
}
inline bool
EDGE_IS_DIRECT(const uint8_t& E)
{
  return ((E >> 1) & 0x3) != 0;
}
inline bool
EDGE_IS_AB(const uint8_t& E)
{
  return ((E >> 1) & 0x1) != 0;
}
inline bool
EDGE_IS_BA(const uint8_t& E)
{
  return ((E >> 2) & 0x1) != 0;
}

inline void
EDGE_SET_EXISTING(uint8_t& E)
{
  (E = (E | 0x1));
}
inline void
EDGE_SET_AB(uint8_t& E)
{
  (E = ((E | 0x2) & 0xFB));
}
inline void
EDGE_SET_BA(uint8_t& E)
{
  (E = ((E | 0x4) & 0xFD));
}
inline void
EDGE_SET_UNDIRECTED(uint8_t& E)
{
  (E = (E & 0xF9));
}
inline void
EDGE_SET_MISSING(uint8_t& E)
{
  (E = (E & 0xFE));
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#ifdef MAKE_EXTERN
extern "C"
{
#endif

  struct __attribute__((packed)) header_attr
  {
    uint64_t nodes = 0;
    uint64_t edges = 0;
    uint16_t nalgs = 0;
    uint8_t dense = 0;
    uint32_t nsupp = 0;
    uint16_t nsupp_str = 0;
    uint16_t nsupp_int = 0;
    uint16_t nsupp_flt = 0;
    uint8_t pagerank_calc = 0;
    uint8_t closeness_calc = 0;
    uint8_t betweenness_calc = 0;
    uint8_t strength_calc = 0;
    uint8_t eigenvector_calc = 0;
    uint8_t katz_calc = 0;

    char version[HEADER_VERSION_SIZE];
    char cmd[HEADER_CMD_SIZE];
  };

  struct __attribute__((packed)) header_node
  {
    char data[HEADER_NODE_SIZE];
  };

  struct __attribute__((packed)) header_alg
  {
    char data[HEADER_ALG_SIZE];
  };

  struct __attribute__((packed)) header_supp
  {
    char data[HEADER_SUPP_SIZE];
  };

  struct __attribute__((packed)) header_centrality
  {
    double data = 0;
  };

  struct __attribute__((packed)) edge_index
  {
    uint32_t i = 0;
    uint32_t j = 0;
  };

  struct __attribute__((packed)) edge_attr
  {
    uint8_t flag = 0;
  };

  struct __attribute__((packed)) edge_score
  {
    double r = 0.0; // rank
    double s = 0.0; // score
  };

  struct __attribute__((packed)) edge_supp_str
  {
    char data[EDGE_SUPP_SIZE];
  };

  struct __attribute__((packed)) edge_supp_int
  {
    int data;
  };

  struct __attribute__((packed)) edge_supp_flt
  {
    float data;
  };

  struct __attribute__((packed)) index_header
  {
    uint8_t dense = 0;
    uint64_t size = 0;
    uint32_t nodes = 0;
    uint64_t edges = 0;
  };

  struct __attribute__((packed)) index_offset
  {
    int64_t data;
  };

  struct __attribute__((packed)) index_index
  {
    uint32_t i = 0;
    uint32_t j = 0;
  };

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  class SeidrFile;
  class SeidrFileHeader;
  class SeidrFileEdge;

  class offset_t
  {
  public:
    uint32_t i = 0;
    uint32_t j = 0;
    int64_t o = 0;
    std::string err;
    bool operator<(const offset_t& rhs) const
    {
      if (i == rhs.i) {
        return j < rhs.j;
      }
      return i < rhs.i;
    }
  };

  class version_t
  {
  public:
    uint32_t major = 0;
    uint32_t minor = 0;
    uint32_t patch = 0;
    version_t() = default;
    version_t(char * V) {
      std::string s;
      uint32_t p = 0;
      for (uint32_t i = 0; i < HEADER_VERSION_SIZE; i++) {
        char c = V[i];
        if (c == '.' || c == '\0') {
          if (p == 0) {
            major = std::stoul(s);
          } else if (p == 1) {
            minor = std::stoul(s);
          } else if (p == 2) {
            patch = std::stoul(s);
          }
          p++;
          s.clear();
          if (c == '\0') {
            break;
          }
          continue;
        }
        s += c;
      }
    }
  };

  class SeidrFile
  {
  public:
    // ctor
    SeidrFile(const char* fp)
      : filepath(fp)
      , bgzfile(nullptr)
      , _closed(0)
      , _opened(0)
    {}
    const char* filepath;
    BGZF* bgzfile;
    void open(const char* mode);
    void close();
    void seek(int64_t offset);
    int64_t tell();
    uint8_t _closed;
    uint8_t _opened;
    void each_edge(std::function<void(SeidrFileEdge&, SeidrFileHeader&)> f,
                   bool include_missing = false);
    void each_edge_exit_early(
      std::function<bool(SeidrFileEdge&, SeidrFileHeader&)> f);
  };

  class SeidrFileHeader
  {
  public:
    header_attr attr;
    version_t version;
    uint16_t get_supp_ind(std::string supp_n);
    bool have_supp(std::string supp_n);
    std::vector<std::string> nodes;
    std::vector<std::string> algs;
    std::vector<std::string> supp;
    void serialize(SeidrFile& file) const;
    void unserialize(SeidrFile& file);
    void print(std::ostream& out);
    void print_centrality(std::ostream& out) const;
    void print(void) { print(std::cout); }
    void print_centrality(void) const { print_centrality(std::cout); }
    std::vector<double> pagerank;
    std::vector<double> closeness;
    std::vector<double> betweenness;
    std::vector<double> strength;
    std::vector<double> eigenvector;
    std::vector<double> katz;
    std::vector<std::string> centrality_names;
    std::vector<double> centrality_data;
    void cmd_from_args(char** argv, int argc);
    void cmd_from_args(const std::vector<std::string>& args,
                       std::string cmd = "");
    void version_from_char(const char* version);
  };

  class SeidrFileEdge
  {
  public:
    edge_index index;
    edge_attr attr;
    std::vector<edge_score> scores;
    std::vector<std::string> supp_str;
    std::vector<int> supp_int;
    std::vector<float> supp_flt;
    void serialize(SeidrFile& file, SeidrFileHeader& header) const;
    void unserialize(SeidrFile& file, SeidrFileHeader& header);
    void print(std::ostream&,
               const SeidrFileHeader& header,
               char end = '\n',
               bool print_supp = true,
               std::string delim = ";",
               std::string sc_delim = ";",
               bool full = false,
               bool supp_tag = true,
               bool no_name = false) const;
    void print(SeidrFileHeader& header) const { print(std::cout, header); }
    void raise_err(std::string what,
                   SeidrFile& file,
                   SeidrFileHeader& header) const;
  };

  class SeidrFileIndex
  {
  public:
    SeidrFileIndex(bool case_insensitive = false, bool strict = false)
      : _case_insensitive(case_insensitive)
      , _strict(strict)
    {}
    void build(std::string& inf);
    void build(SeidrFile& inf);
    void serialize(SeidrFile& file);
    void unserialize(SeidrFile& file, SeidrFileHeader& head);
    uint64_t coord_to_inx(uint32_t i, uint32_t j)
    {
      return (((i + 1) * i) / 2) - ((i + 1) - (j + 1));
    }
    std::pair<bool, uint32_t> find(const std::string& x);
    std::vector<offset_t> get_offset_node(std::string& node);
    offset_t get_offset_pair(std::string& lhs, std::string& rhs);
    std::set<offset_t> get_offset_nodelist(std::vector<std::string>& nodelist);
    //
    std::vector<int64_t> data;
    uint8_t dense;
    uint32_t nodes;
    uint64_t edges;
    std::map<std::string, uint32_t> gmap;

  private:
    bool _case_insensitive;
    bool _strict;
  };

  struct MiniEdge
  {
    uint64_t i = 0;
    uint64_t j = 0;
    double s = 0;
  };

  bool me_score_sort(MiniEdge a, MiniEdge b);
  bool me_rank_sort(MiniEdge a, MiniEdge b);
  bool me_score_sort_abs(MiniEdge a, MiniEdge b);
  bool me_score_sort_rev(MiniEdge a, MiniEdge b);

  bool sfe_score_sort(SeidrFileEdge a, SeidrFileEdge b);
  bool sfe_rank_sort(SeidrFileEdge a, SeidrFileEdge b);

  class sfe_rank_sort_class
  {
  public:
    bool operator()(const SeidrFileEdge& lhs, const SeidrFileEdge& rhs) const
    {
      return sfe_rank_sort(lhs, rhs);
    }
  };

  class sfe_score_sort_class
  {
  public:
    bool operator()(const SeidrFileEdge& lhs, const SeidrFileEdge& rhs) const
    {
      return sfe_score_sort(lhs, rhs);
    }
  };

  std::vector<MiniEdge> read_network_minimal(
    SeidrFileHeader& h,
    SeidrFile& rf,
    double threshold = -std::numeric_limits<double>::infinity(),
    uint32_t tpos = 0,
    bool trank = false);

  arma::mat read_network_arma(
    SeidrFileHeader& h,
    SeidrFile& rf,
    double threshold = -std::numeric_limits<double>::infinity(),
    uint32_t tpos = 0,
    bool trank = false);

  std::vector<SeidrFileEdge> read_network(
    SeidrFileHeader& h,
    SeidrFile& rf,
    double threshold = -std::numeric_limits<double>::infinity(),
    uint32_t tpos = 0,
    bool trank = false);

  class SeidrFileIndexNodeNotFound : public std::exception
  {
  public:
    SeidrFileIndexNodeNotFound(std::string msg)
      : _msg(msg)
    {}

    virtual const char* what() const throw() { return _msg.c_str(); }

  private:
    std::string _msg;
  };

#ifdef MAKE_EXTERN
}
#endif