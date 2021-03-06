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
#include <BSlogger.hpp>
#include <Serialize.h>
#include <common.h>
#include <fs.h>
#include <import.h>
#include <parallel_control.h>
// External
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Parallel includes
#if defined(SEIDR_PSTL)
#include <pstl/algorithm>
#include <pstl/execution>
#include <tbb/task_scheduler_init.h>
#else
#include <algorithm>
#endif

#include <armadillo>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <map>
#include <memory>

namespace po = boost::program_options;
using boost::lexical_cast;
using boost::numeric_cast;
using stringmap = std::map<std::string, uint32_t>;
using inx_t = std::pair<uint32_t, uint32_t>;

seidr_score_t
parse_score_field(const std::string& field, const uint64_t& l)
{
  try {
    return std::stod(field);
  } catch (std::exception& e) {
    std::stringstream ss;
    ss << "Exception in line: [" << l << "], tried to parse numeric: [" << field
       << "]. Check file input format. Exception: " << e.what();
    throw std::runtime_error(ss.str());
  }
}

read_logger&
read_logger::operator++()
{
  _i++;
  if (_i % _n == 0) {
    _l(LOG_INFO) << "Read " << _i << " edges\n";
  }
  return *this;
}

read_logger // NOLINT: Overlaps with readability-const-return-type
read_logger::operator++(int)
{
  read_logger tmp = *this;
  ++*this;
  return tmp;
}

void
lt_map::reserve(uint64_t n)
{
  //_ev = std::vector<edge>(n);
  _shadow = std::vector<shadow_t>(n);
  // Make sure default ctor is called
  for (auto& i : _shadow) {
    i = shadow_t();
  }
}

void
lt_map::insert(inx_t& inx, reduced_edge& ee, bool& rev, bool& drop)
{
  // Always work in lower triangle
  if (inx.first < inx.second) {
    uint32_t tmp = inx.first;
    inx.first = inx.second;
    inx.second = tmp;
    // Keep track of edge direction
    ee.d = 2;
  } else if (ee.d != 3) {
    ee.d = 1;
  } else {
    ee.d = 0;
  }
  edge e;
  e.i = inx.first;
  e.j = inx.second;
  e.d = ee.d;
  e.w = ee.w;

  uint64_t shadow_offset = _hasher(e);
  // Check if edge already exists
  shadow_t* ptr = &_shadow[shadow_offset];
  // If not, insert edge and do nothing else
  if (ptr->found == 0) {
    if (!drop || (drop && (!almost_equal(e.w, 0)))) {
      _ev.push_back(e);
      ptr->offset = _ev.size() - 1;
      ptr->found = 1;
    }
  }
  // If there is a symmetric edge present, check if current edge has a better
  // score and insert if it has.
  else if (!drop || (drop && !almost_equal(e.w, 0))) {
    // In reverse mode, bigger is better
    edge* ev_ptr = &_ev[ptr->offset];
    if (rev) {
      if (e.w > ev_ptr->w) {
        ev_ptr->w = e.w;
        ev_ptr->d = e.d;
      } else if (almost_equal(e.w, ev_ptr->w)) {
        // If A->B == B<-A, set to undirected.
        ev_ptr->d = 0;
      }
    }
    // In default mode, lower is better
    else {
      if (e.w < ev_ptr->w) {
        ev_ptr->w = e.w;
        ev_ptr->d = e.d;
      } else if (almost_equal(e.w, ev_ptr->w)) {
        // If A->B == B<-A, set to undirected.
        ev_ptr->d = 0;
      }
    }
  }
}

// Convert lower triangular map to vector of edges
std::vector<edge>&
lt_map::to_vec()
{
  return _ev;
}

/**
 * A comparison to rank scores in ascending or descending order
 */
bool
ascending(const edge& a, const edge& b)
{
  return (a.w < b.w);
}
bool
descending(const edge& a, const edge& b)
{
  return (a.w > b.w);
}
bool
abs_ascending(const edge& a, const edge& b)
{
  return (fabs(a.w) < fabs(b.w));
}
bool
abs_descending(const edge& a, const edge& b)
{
  return (fabs(a.w) > fabs(b.w));
}
bool
lt(const edge& a, const edge& b)
{
  if (a.i < b.i) {
    return true;
  }
  if (a.i == b.i) {
    return a.j < b.j;
  }
  return false;
}
/**
 * A function to replace the scores in a vector with the
 * ranks they assume in the entire set.
 *
 * @param a A vector iterator pointing to the start of the vector.
 * @param b A vector iterator pointing to the end of the vector.
 * @param reverse Boolean indicator to check if the vector should
 *                be sorted in descending order.
 */
void
rank_vector(std::vector<edge>& ev, bool reverse, bool absolute)
{
  // Sorting
  logger log(std::cerr, "rank_vector");

  // Do nothing on empty dataset
  if (ev.empty()) {
    return;
  }

  if (reverse) {
    if (absolute) {
      SORTWCOMP(ev.begin(), ev.end(), abs_descending);
    } else {
      SORTWCOMP(ev.begin(), ev.end(), descending);
    }
  } else {
    if (absolute) {
      SORTWCOMP(ev.begin(), ev.end(), abs_ascending);
    } else {
      SORTWCOMP(ev.begin(), ev.end(), ascending);
    }
  }
  auto it = ev.begin();
  uint64_t pos = 0;
  seidr_score_t prev = it->w;
  uint64_t start = 0;
  seidr_score_t rank = 0;
  while (it != ev.end()) {
    it++;
    pos++;
    if (it == ev.end() || it->w != prev) {
      rank = (lexical_cast<seidr_score_t>(pos) + 1 +
              lexical_cast<seidr_score_t>(start)) /
             2;
      for (uint64_t i = start; i < pos; i++) {
        // edge e = ev[i];
        // e.r = rank;
        ev[i].r = rank;
      }
      if (it != ev.end()) {
        start = pos;
        prev = it->w;
      }
    }
  }
  SORTWCOMP(ev.begin(), ev.end(), lt);
}

/**
 * Creating a binary ranked file from an edge list input.
 *
 * @param el_file The input edge list.
 * @param gene_file The gene names.
 * @return int 0 if the function succeeded, an error code otherwise.
 */
int
import(const std::vector<std::string>& args)
{

  logger log(std::cerr, "import");

  try {

    read_logger pr(log, 100000000); // NOLINT

    std::vector<std::string> genes;

    // Variables used by the function
    seidr_import_param_t param;

    po::options_description umbrella(
      "Convert various text based network representations to "
      "SeidrFiles");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message");

    po::options_description ropt("Ranking Options");
    ropt.add_options()("absolute,A",
                       po::bool_switch(&param.absolute)->default_value(false),
                       "Rank on absolute of scores")(
      "reverse,r",
      po::bool_switch(&param.rev)->default_value(false),
      "Rank scores in descending order (highest first)")(
      "drop-zero,z",
      po::bool_switch(&param.drop_zero)->default_value(false),
      "Drop edges with a weight of zero from the network")(
      "undirected,u",
      po::bool_switch(&param.force_undirected)->default_value(false),
      "Force all edges to be interpreted as undirected");

    po::options_description ompopt("OpenMP Options");
    ompopt.add_options()(
      "threads,O",
      po::value<int>(&param.nthreads)->default_value(GET_MAX_PSTL_THREADS()),
      "Number of OpenMP threads for parallel sorting");

    po::options_description req("Required");
    req.add_options()("infile,i",
                      po::value<std::string>(&param.el_file)->required(),
                      "Input file name ['-' for stdin]")(
      "format,F",
      po::value<std::string>(&param.format)->required(),
      "The input file format [el, lm, m, ara]")(
      "genes,g",
      po::value<std::string>(&param.gene_file)->required(),
      "Gene file for input")(
      "outfile,o",
      po::value<std::string>(&param.out_file)->required(),
      "Output file name")("name,n",
                          po::value<std::string>(&param.name)->required(),
                          "Import name (algorithm name)");

    umbrella.add(req).add(ropt).add(ompopt).add(opt);

    po::variables_map vm;
    po::store(po::command_line_parser(args).options(umbrella).run(), vm);

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);

    set_pstl_threads(param.nthreads);

    if (param.format == "el") {
      param.is_edge_list = true;
    } else if (param.format == "lm") {
      param.is_lm = true;
    } else if (param.format == "m") {
      param.is_full_matrix = true;
    } else if (param.format == "ara") {
      param.is_aracne = true;
    } else {
      throw std::runtime_error("Unknown input format: " + param.format);
    }

    if (param.el_file != "-") {
      param.el_file = to_absolute(param.el_file);
      assert_exists(param.el_file);
      assert_can_read(param.el_file);
    }

    param.gene_file = to_absolute(param.gene_file);
    assert_exists(param.gene_file);
    assert_can_read(param.gene_file);

    param.out_file = to_absolute(param.out_file);
    if (!param.force) {
      assert_no_overwrite(param.out_file);
    }
    assert_dir_is_writeable(dirname(param.out_file));
    assert_arg_constraint<std::string>({ "el", "lm", "m", "ara" },
                                       param.format);

    genes = read_genes(param.gene_file, '\n', '\t');
    stringmap g_map;

    size_t i = 0;
    for (std::string& s : genes) {
      g_map.insert(std::pair<std::string, size_t>(s, i++));
    }

    log(LOG_INFO) << "Read " << genes.size() << " gene names\n";

    log(LOG_INFO) << "Allocating initial matrix\n";
    lt_map X;

    X.reserve((genes.size() * (genes.size() - 1)) / 2);

    log(LOG_INFO) << "Populating matrix\n";
    std::shared_ptr<std::istream> ifs;

    if (param.el_file == "-") {
      ifs.reset(&(std::cin), no_delete);
    } else {
      ifs = std::shared_ptr<std::istream>(
        new std::ifstream(param.el_file.c_str(), std::ifstream::in));
    }

    uint64_t lno = 0;

    if (param.is_edge_list) {
      //    std::istream ifs = *in_stream;
      if (!(*ifs).good()) {
        throw std::runtime_error(
          "Error opening edge list. Check permissions/spelling.");
      }
      seidr_score_t v = 0; // value
      std::string gi;
      std::string gj;

      for (std::string row; std::getline((*ifs), row, '\n');) {
        lno++;
        std::istringstream fl(row);
        size_t x = 0;
        for (std::string field; std::getline(fl, field, '\t');) {
          if (x == 0) {
            gi = field;
          }
          if (x == 1) {
            gj = field;
          }
          if (x == 2) {
            v = parse_score_field(field, lno);
          }
          x++;
        }
        reduced_edge e;
        inx_t inx;
        try {
          inx.first = g_map.at(gi);
          inx.second = g_map.at(gj);
        } catch (std::exception& e) {
          std::stringstream ss;
          ss << "Cannot find key pair: [" << gi << ", " << gj
             << "] in gene map. These should be gene names in an edge-list "
                "formatted "
             << "file.";
          throw std::runtime_error(ss.str());
        }
        e.w = v;
        X.insert(inx, e, param.rev, param.drop_zero);
        pr++;
      }
      if ((*ifs).bad()) {
        throw std::runtime_error("Error reading edge list.");
      }
      //(*ifs).close();
    } else if (param.is_lm) {
      // Lower triangular matrix, no diagonal
      //      std::istream ifs = *in_stream;
      if (!(*ifs).good()) {
        throw std::runtime_error(
          "Error opening edge list. Check permissions/spelling.");
      }
      for (uint32_t i = 1; i < genes.size(); i++) {
        std::string line;
        std::getline((*ifs), line);
        lno++;
        std::stringstream ss(line);
        for (uint32_t j = 0; j < i; j++) {
          std::string col;
          ss >> col;
          seidr_score_t x = parse_score_field(col, lno);
          reduced_edge e;
          inx_t inx;
          inx.first = i;
          inx.second = j;
          e.w = x;
          X.insert(inx, e, param.rev, param.drop_zero);
          pr++;
        }
      }
    } else if (param.is_full_matrix) {
      if (!(*ifs).good()) {
        throw std::runtime_error(
          "Error opening edge list. Check permissions/spelling.");
      }
      for (uint32_t i = 0; i < genes.size(); i++) {
        std::string line;
        std::getline((*ifs), line);
        lno++;
        std::stringstream ss(line);
        for (uint32_t j = 0; j < genes.size(); j++) {
          std::string col;
          ss >> col;
          if (i == j) {
            continue;
          }
          seidr_score_t x = parse_score_field(col, lno);
          reduced_edge e;
          inx_t inx;
          inx.first = i;
          inx.second = j;
          e.w = x;
          X.insert(inx, e, param.rev, param.drop_zero);
          pr++;
        }
      }
    } else if (param.is_aracne) {
      // ARACNE code
      //      std::istream ifs = *in_stream;
      if (!(*ifs).good()) {
        throw std::runtime_error(
          "Error opening edge list. Check permissions/spelling.");
      }
      std::string line;
      while (std::getline((*ifs), line)) {
        lno++;
        if (line.at(0) == '>') {
          continue; // comments
        }
        std::stringstream ss(line);
        std::string gi;
        ss >> gi;
        uint32_t i = g_map.at(gi);
        uint32_t tmp = 1;
        uint32_t j = 0;
        seidr_score_t x = 0;
        std::string col;
        std::string target;
        while (ss.good()) {
          ss >> col;
          if (tmp % 2 == 1) {
            target = col;
          } else if (tmp % 2 == 0) {
            j = g_map.at(target);
            x = parse_score_field(col, lno);
            reduced_edge e;
            inx_t inx;
            inx.first = i;
            inx.second = j;
            e.w = x;
            X.insert(inx, e, param.rev, param.drop_zero);
            pr++;
          }
          tmp++;
        }
      }
    } else {
      // Shouldn't happen...
      throw std::runtime_error("Unknown file format.");
    }

    // log(LOG_INFO) << "Read " << pr << " edges\n";

    log(LOG_INFO) << "Swapping min/max edges to lower triangle\n";
    std::vector<edge>& vs = X.to_vec();

    log(LOG_INFO) << "Computing ranks\n";
    rank_vector(vs, param.rev, param.absolute);

    log(LOG_INFO) << "Writing data: " << genes.size() << " nodes, " << vs.size()
                  << " edges\n";

    uint64_t cnt = 0;
    uint8_t dense = 0;
    double nn = genes.size();
    double ne = vs.size();
    if (ne >
        ((nn * (nn - 1) / 2) * 0.66)) { // NOLINT(readability-magic-numbers)
      dense = 1;
    }

    SeidrFile ostr(param.out_file.c_str());
    ostr.open("w");

    SeidrFileHeader h;
    h.attr.nodes = genes.size();
    h.attr.edges = vs.size();
    h.attr.nalgs = 1;
    h.attr.dense = dense;
    h.version_from_char(_XSTR(VERSION));
    h.cmd_from_args(args);

    for (const auto& it : genes) {
      h.nodes.push_back(it);
    }

    std::vector<std::string> method = { param.name };
    h.algs = method;
    h.serialize(ostr);

    // Main data. di and dj track the lower triangular index
    // of the current edge and empty edges will be created
    // in case they are missing in the rank vector. In sparse
    // mode these are ignored
    uint32_t di = 1;
    uint32_t dj = 0;
    uint64_t end =
      dense != 0 ? (genes.size() * (genes.size() - 1)) / 2 : vs.size();
    for (size_t i = 0; i < end; i++) {
      SeidrFileEdge e;

      // Create an empty entry for missing edges in dense
      // storage mode
      if (dense != 0 && i < vs.size()) {
        while (vs[i].i != di || vs[i].j != dj) {
          SeidrFileEdge tmp;
          EDGE_SET_MISSING(tmp.attr.flag);
          tmp.serialize(ostr, h);
          if (dj == di - 1) {
            di++;
            dj = 0;
          } else {
            dj++;
          }
        }
        // If i,j of the lower triangular index are euqal
        // to the current node, increment them to the next
        // node in the index for the next iteration
        if (dj == di - 1) {
          di++;
          dj = 0;
        } else {
          dj++;
        }
      }

      // Only relevant in sparse mode
      if (dense == 0) {
        e.index.i = vs[i].i;
        e.index.j = vs[i].j;
      }
      if (i < vs.size()) {
        edge_score s;
        s.s = vs[i].w;
        s.r = vs[i].r;
        e.scores.push_back(s);

        EDGE_SET_EXISTING(e.attr.flag);

        if (param.force_undirected) {
          EDGE_SET_UNDIRECTED(e.attr.flag);
        } else {
          switch (vs[i].d) {
            case 0:
              EDGE_SET_UNDIRECTED(e.attr.flag);
              break;
            case 1:
              EDGE_SET_AB(e.attr.flag);
              break;
            case 2:
              EDGE_SET_BA(e.attr.flag);
              break;
            default:
              throw std::runtime_error("Unknown edge direction\n");
          }
        }
        e.serialize(ostr, h);
        cnt++;
        if (cnt % PRINTING_MOD == 0) {
          log(LOG_INFO) << cnt << " edges written ("
                        << (numeric_cast<double>(cnt) /
                            numeric_cast<double>(vs.size())) *
                             100.0
                        << "%)\n";
        }
      } else {
        SeidrFileEdge tmp;
        EDGE_SET_MISSING(tmp.attr.flag);
        tmp.serialize(ostr, h);
      }
    }
    ostr.close();
  } catch (const po::error& e) {
    log(LOG_ERR) << "[Argument Error]: " << e.what() << '\n';
    return 1;
  } catch (const std::runtime_error& e) {
    log(LOG_ERR) << "[Runtime Error]: " << e.what() << '\n';
    return 1;
  } catch (const std::exception& e) {
    log(LOG_ERR) << "[Generic Error]: " << e.what() << '\n';
    return 1;
  }

  log.time_since_start();
  return 0;
}
