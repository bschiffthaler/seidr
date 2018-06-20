#pragma once

#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <htslib/bgzf.h>
#include <map>
#include <set>
#include <limits>
#include <armadillo>

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
#define EDGE_EXISTS(E) (E & 0x1)
#define EDGE_IS_DIRECT(E) ( (E>>1) & 0x3 )
#define EDGE_IS_AB(E) ( (E>>1) & 0x1 )
#define EDGE_IS_BA(E) ( (E>>2) & 0x1 )
#define EDGE_SET_EXISTING(E) (E = (E | 0x1))
#define EDGE_SET_AB(E) (E = ((E | 0x2) & 0xFB))
#define EDGE_SET_BA(E) (E = ((E | 0x4) & 0xFD))
#define EDGE_SET_UNDIRECTED(E) (E = (E & 0xF9))
#define EDGE_SET_MISSING(E) (E = (E & 0xFE))

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

struct __attribute__ ((packed)) header_attr {
  uint64_t nodes = 0;
  uint64_t edges = 0;
  uint16_t nalgs = 0;
  uint8_t  dense = 0;
  uint32_t nsupp = 0;
  uint16_t nsupp_str = 0;
  uint16_t nsupp_int = 0;
  uint16_t nsupp_flt = 0;
  uint8_t  pagerank_calc = 0;
  uint8_t  closeness_calc = 0;
  uint8_t  betweenness_calc = 0;
  uint8_t  strength_calc = 0;
  uint8_t  eigenvector_calc = 0;
  uint8_t  katz_calc = 0;

  char version[HEADER_VERSION_SIZE];
  char         cmd[HEADER_CMD_SIZE];
};

struct __attribute__ ((packed)) header_node {
  char data[HEADER_NODE_SIZE];
};

struct __attribute__ ((packed)) header_alg {
  char data[HEADER_ALG_SIZE];
};

struct __attribute__ ((packed)) header_supp {
  char data[HEADER_SUPP_SIZE];
};

struct __attribute__ ((packed)) header_centrality {
  double data = 0;;
};

struct __attribute__ ((packed)) edge_index {
  uint32_t i = 0;
  uint32_t j = 0;
};

struct __attribute__ ((packed)) edge_attr {
  uint8_t flag = 0;
};

struct __attribute__ ((packed)) edge_score {
  double r = 0.0; //rank
  double s = 0.0; //score
};

struct __attribute__ ((packed)) edge_supp_str {
  char data[EDGE_SUPP_SIZE];
};

struct __attribute__ ((packed)) edge_supp_int {
  int data;
};

struct __attribute__ ((packed)) edge_supp_flt {
  float data;
};

struct __attribute__ ((packed)) index_header {
  uint8_t dense = 0;
  uint64_t size = 0;
  uint32_t nodes = 0;
  uint64_t edges = 0;
};

struct __attribute__ ((packed)) index_offset {
  int64_t data;
};

struct __attribute__ ((packed)) index_index {
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
  bool operator<(const offset_t& rhs) const {
    if (i == rhs.i)
      return j < rhs.j;
    else
      return i < rhs.i;
  }
};

class SeidrFile {
public:
//ctor
  SeidrFile(const char * fp) : filepath(fp), _closed(0), _opened(0) {}
  const char * filepath;
  BGZF* bgzfile;
  void open(const char * mode);
  void close();
  void seek(int64_t offset);
  uint8_t _closed;
  uint8_t _opened;
  void each_edge(std::function<void(SeidrFileEdge&, SeidrFileHeader&)> f);
};

class SeidrFileHeader {
public:
  header_attr attr;
  uint16_t get_supp_ind(std::string supp_n);
  bool have_supp(std::string supp_n);
  std::vector<std::string> nodes;
  std::vector<std::string> algs;
  std::vector<std::string> supp;
  void serialize(SeidrFile& file);
  void unserialize(SeidrFile& file);
  void print(std::ostream& out);
  void print_centrality(std::ostream& out);
  void print(void) {print(std::cout);}
  void print_centrality(void) {print_centrality(std::cout);}
  std::vector<double> pagerank;
  std::vector<double> closeness;
  std::vector<double> betweenness;
  std::vector<double> strength;
  std::vector<double> eigenvector;
  std::vector<double> katz;
  void cmd_from_args(char ** argv, int argc);
  void version_from_char(const char * version);
};

class SeidrFileEdge {
public:
  edge_index index;
  edge_attr   attr;
  std::vector<edge_score> scores;
  std::vector<std::string> supp_str;
  std::vector<int> supp_int;
  std::vector<float> supp_flt;
  void serialize(SeidrFile& file, SeidrFileHeader& header);
  void unserialize(SeidrFile& file, SeidrFileHeader& header);
  void print(std::ostream&, SeidrFileHeader& header,  char end = '\n',
             bool print_supp = true, std::string delim = ";",
             std::string sc_delim = ";", bool full = false,
             bool supp_tag = true, bool no_name = false) const;
  void print(SeidrFileHeader& header) const {print(std::cout, header);}
  void raise_err(std::string what, SeidrFile& file, SeidrFileHeader& header);
};

class SeidrFileIndex {
public:
  void build(std::string& inf);
  void build(SeidrFile& inf);
  void serialize(SeidrFile& file);
  void unserialize(SeidrFile& file, SeidrFileHeader& head);
  uint64_t coord_to_inx(uint32_t i, uint32_t j) {
    return (((i + 1) * i) / 2) - ((i + 1) - (j + 1));
  }
  uint32_t find(std::string& x);
  std::vector<offset_t> get_offset_node(std::string& node);
  offset_t get_offset_pair(std::string& lhs, std::string& rhs);
  std::set<offset_t> get_offset_nodelist(std::vector<std::string>& nodelist);
  //
  std::vector<int64_t> data;
  uint8_t dense;
  uint32_t nodes;
  uint64_t edges;
  std::map<std::string, uint32_t> gmap;
};

struct MiniEdge {
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

class sfe_rank_sort_class {
public:
    bool operator()(const SeidrFileEdge& lhs, const SeidrFileEdge& rhs) const { 
        return sfe_rank_sort(lhs, rhs); 
    }
};

class sfe_score_sort_class {
public:
    bool operator()(const SeidrFileEdge& lhs, const SeidrFileEdge& rhs) const { 
        return sfe_score_sort(lhs, rhs); 
    }
};

std::vector<MiniEdge>
read_network_minimal(SeidrFileHeader& h,
                     SeidrFile& rf,
                     double threshold = -std::numeric_limits<double>::infinity(),
                     uint32_t tpos = 0,
                     bool trank = false);

arma::mat
read_network_arma(SeidrFileHeader& h,
                  SeidrFile& rf,
                  double threshold = -std::numeric_limits<double>::infinity(),
                  uint32_t tpos = 0,
                  bool trank = false);

std::vector<SeidrFileEdge>
read_network(SeidrFileHeader& h,
             SeidrFile& rf,
             double threshold = -std::numeric_limits<double>::infinity(),
             uint32_t tpos = 0,
             bool trank = false);
