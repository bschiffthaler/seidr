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

// Seir
#include <BSlogger.hpp>
#include <Serialize.h>
#include <aggregate.h>
#include <common.h>
#include <fs.h>
#include <parallel_control.h>
// Parallel includes
#if defined(SEIDR_PSTL)
#include <pstl/algorithm>
#include <pstl/execution>
#else
#include <algorithm>
#endif
// External
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <cerrno>
#include <cmath>
#include <forward_list>
#include <fstream>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace po = boost::program_options;
using boost::lexical_cast;
using boost::numeric_cast;
namespace fs = boost::filesystem;

struct aggr_rank_t
{
  double score = 0;
  double rank = 0;
  uint64_t index = 0;
};

bool
ev_score_sort(aggr_rank_t a, aggr_rank_t b)
{
  return (a.score < b.score);
}
bool
ev_index_sort(aggr_rank_t a, aggr_rank_t b)
{
  return (a.index < b.index);
}

uint8_t
majority_dir_vote(std::vector<uint8_t>& flags)
{
  uint8_t ret = 0;
  uint16_t ab = 0;
  uint16_t ba = 0;
  for (auto& flag : flags) {
    // Ignore non-directional data
    if (!EDGE_IS_DIRECT(flag)) {
      continue;
    }
    if (EDGE_IS_AB(flag)) {
      ab++;
    } else {
      ba++;
    }
  }
  if (ab > (2 * ba)) {
    EDGE_SET_AB(ret);
  } else if (ba > (2 * ab)) {
    EDGE_SET_BA(ret);
  }
  return ret;
}

void
rank_vector(std::vector<aggr_rank_t>& ev)
{
  logger log(std::cerr, "rank_vector");
  SORTWCOMP(ev.begin(), ev.end(), ev_score_sort);
  log(LOG_DEBUG) << "Computing ranks\n";
  auto it = ev.begin();
  uint64_t pos = 0;
  double prev = it->score;
  uint64_t start = 0;
  double rank;
  while (it != ev.end()) {
    it++;
    pos++;
    if (it == ev.end() || !almost_equal(it->score, prev)) {
      rank = (lexical_cast<double>(pos) + 1 + lexical_cast<double>(start)) / 2;
      for (size_t i = start; i < pos; i++) {
        ev[i].rank = rank;
      }
      if (it != ev.end()) {
        start = pos;
        prev = it->score;
      }
    }
  }
  SORTWCOMP(ev.begin(), ev.end(), ev_index_sort);
}

void
aggr_top1(std::vector<SeidrFileHeader>&
            header_vec, // NOLINT(clang-diagnostic-unused-parameter)
          std::vector<SeidrFileEdge>& v,
          SeidrFileEdge& result,
          std::vector<uint8_t>& existing,
          std::vector<uint8_t>& flags,
          double& rank,
          uint64_t& ex_sum,
          bool& keep_di)
{
  uint16_t index = 0;
  for (uint16_t index_b = 0;
       index_b < numeric_cast<uint16_t>(result.scores.size());
       index_b++) {
    if (existing[index_b] != 0 && result.scores[index_b].r < rank) {
      rank = result.scores[index_b].r;
      index = index_b;
    }
  }

  edge_score es;
  es.r = rank;
  es.s = 0;
  result.scores.push_back(es);

  result.supp_int.push_back(numeric_cast<uint16_t>(index));

  if (keep_di) {
    for (auto& flag : flags) {
      result.supp_int.push_back(numeric_cast<uint16_t>(flag));
    }
  }

  if (ex_sum > 0) {
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_DIRECT(v[index].attr.flag)) {
      if (EDGE_IS_AB(v[index].attr.flag)) {
        EDGE_SET_AB(result.attr.flag);
      } else {
        EDGE_SET_BA(result.attr.flag);
      }
    }
  }
}

void
aggr_top2(std::vector<SeidrFileHeader>&
            header_vec, // NOLINT(clang-diagnostic-unused-parameter)
          std::vector<SeidrFileEdge>& v,
          SeidrFileEdge& result,
          std::vector<uint8_t>& existing,
          std::vector<uint8_t>& flags,
          double& rank,
          uint64_t& ex_sum,
          bool& keep_di)
{

  uint16_t index1 = 0;
  uint16_t index2 = 0;
  double r1 = rank;
  double r2 = rank;
  for (uint16_t index_b = 0;
       index_b < numeric_cast<uint16_t>(result.scores.size());
       index_b++) {
    if (existing[index_b] != 0 && result.scores[index_b].r < r1) {
      r1 = result.scores[index_b].r;
      index1 = index_b;
    } else if (existing[index_b] != 0 && result.scores[index_b].r < r2) {
      r2 = result.scores[index_b].r;
      index2 = index_b;
    }
  }

  edge_score es;
  es.r = (r1 + r2) / 2;
  es.s = 0;
  result.scores.push_back(es);

  result.supp_int.push_back(numeric_cast<uint16_t>(index1));
  result.supp_int.push_back(numeric_cast<uint16_t>(index2));

  if (keep_di) {
    for (auto& flag : flags) {
      result.supp_int.push_back(numeric_cast<uint16_t>(flag));
    }
  }

  if (ex_sum > 0) {
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_DIRECT(v[index1].attr.flag) &&
        EDGE_IS_DIRECT(v[index2].attr.flag)) {
      if (EDGE_IS_AB(v[index1].attr.flag) && EDGE_IS_AB(v[index2].attr.flag)) {
        EDGE_SET_AB(result.attr.flag);
      } else if (EDGE_IS_BA(v[index1].attr.flag) &&
                 EDGE_IS_BA(v[index2].attr.flag)) {
        EDGE_SET_BA(result.attr.flag);
      }
    }
  }
}

void
aggr_borda(std::vector<SeidrFileHeader>&
             header_vec, // NOLINT(clang-diagnostic-unused-parameter)
           std::vector<SeidrFileEdge>& v,
           SeidrFileEdge& result,
           std::vector<uint8_t>& existing,
           std::vector<uint8_t>& flags,
           double& rank, // NOLINT(clang-diagnostic-unused-parameter)
           uint64_t& ex_sum,
           bool& keep_di)
{
  uint16_t index = 0;
  double sum = 0;
  double cnt = 0;
  for (uint16_t index_b = 0;
       index_b < numeric_cast<uint16_t>(result.scores.size());
       index_b++) {
    if (existing[index_b] != 0) {
      sum += result.scores[index_b].r;
      cnt++;
    }
  }

  edge_score es;
  es.r = sum / cnt;
  es.s = 0;
  result.scores.push_back(es);

  if (keep_di) {
    for (auto& flag : flags) {
      result.supp_int.push_back(numeric_cast<uint16_t>(flag));
    }
  }

  if (ex_sum > 0) {
    uint8_t dir = majority_dir_vote(flags);
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_AB(dir)) {
      EDGE_SET_AB(result.attr.flag);
    } else if (EDGE_IS_BA(dir)) {
      EDGE_SET_BA(result.attr.flag);
    }
  }
}

void
aggr_irp(
  std::vector<SeidrFileHeader>& header_vec,
  std::vector<SeidrFileEdge>& v, // NOLINT(clang-diagnostic-unused-parameter)
  SeidrFileEdge& result,
  std::vector<uint8_t>& existing,
  std::vector<uint8_t>& flags, // NOLINT(clang-diagnostic-unused-parameter)
  double& rank, // NOLINT(clang-diagnostic-unused-parameter)
  uint64_t& ex_sum, // NOLINT(clang-diagnostic-unused-parameter)
  bool& keep_di)
{
  // uint16_t index = 0;
  double r;
  double nodes = header_vec[0].attr.nodes;
  double nmax = log10((nodes * (nodes - 1)) / 2);
  uint16_t cnt = 0;
  if (existing[0] != 0) {
    // As LHS and RHS are logs multiplication of both
    // becomes a sum
    r = 1 + log10(result.scores[0].r);
    cnt++;
  } else {
    r = 1 + nmax;
  }
  for (uint16_t index_b = 1;
       index_b < numeric_cast<uint16_t>(result.scores.size());
       index_b++) {
    if (existing[index_b] != 0) {
      r += (1 + log10(result.scores[index_b].r));
      cnt++;
    } else {
      r += (1 + nmax);
    }
  }

  edge_score es;
  es.r = r;
  es.s = 0;
  result.scores.push_back(es);

  if (keep_di) {
    for (auto& flag : flags) {
      result.supp_int.push_back(numeric_cast<uint16_t>(flag));
    }
  }

  if (cnt > 0) {
    uint8_t dir = majority_dir_vote(flags);
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_AB(dir)) {
      EDGE_SET_AB(result.attr.flag);
    } else if (EDGE_IS_BA(dir)) {
      EDGE_SET_BA(result.attr.flag);
    }
  }
}

SeidrFileEdge
calc_score(std::vector<SeidrFileHeader>& header_vec,
           std::vector<SeidrFileEdge>& v,
           std::vector<uint8_t>& get_next,
           uint32_t i,
           uint64_t j,
           const std::function<void(std::vector<SeidrFileHeader>&,
                                    std::vector<SeidrFileEdge>&,
                                    SeidrFileEdge&,
                                    std::vector<uint8_t>&,
                                    std::vector<uint8_t>&,
                                    double&,
                                    uint64_t&,
                                    bool&)>& aggr_fun,
           bool& keep_di)
{
  SeidrFileEdge result;
  result.index.i = i;
  result.index.j = j;
  std::vector<uint8_t> existing;
  std::vector<uint8_t> flags;
  uint64_t ex_sum = 0;
  double rank = std::numeric_limits<double>::infinity();
  uint32_t index_a = 0;
  edge_score dummy;
  dummy.r = std::numeric_limits<double>::quiet_NaN();
  dummy.s = std::numeric_limits<double>::quiet_NaN();
  for (auto& e : v) {
    flags.push_back(e.attr.flag);
    if (e.index.i == i && e.index.j == j) {
      if (EDGE_EXISTS(e.attr.flag)) {
        for (auto& s : e.scores) {
          existing.push_back(1);
          ex_sum += 1;
          result.scores.push_back(s);
        }
        for (auto& s : e.supp_str) {
          result.supp_str.push_back(s);
        }
        for (auto& s : e.supp_int) {
          result.supp_int.push_back(s);
        }
        for (auto& s : e.supp_flt) {
          result.supp_flt.push_back(s);
        }
      } else {
        existing.push_back(0);
        result.scores.push_back(dummy);
        for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_str; a++) {
          result.supp_str.emplace_back("");
        }
        for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_int; a++) {
          result.supp_int.push_back(0);
        }
        for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_flt; a++) {
          result.supp_flt.push_back(0);
        }
      }
      get_next[index_a] = 1;
    } else
    // Storage is sparse and rank file is ahead of current edge
    // We fill the edge with dummy data
    {
      for (uint16_t a = 0; a < header_vec[index_a].attr.nalgs; a++) {
        result.scores.push_back(dummy);
        existing.push_back(0);
      }
      for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_str; a++) {
        result.supp_str.emplace_back("");
      }
      for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_int; a++) {
        result.supp_int.push_back(0);
      }
      for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_flt; a++) {
        result.supp_flt.push_back(0);
      }
      // Do not read next edge from this file
      get_next[index_a] = 0;
    }
    index_a++;
  }

  // Call function pointer
  aggr_fun(header_vec, v, result, existing, flags, rank, ex_sum, keep_di);

  return result;
}

int
aggregate(const std::vector<std::string>& args)
{

  logger log(std::cerr, "aggregate");

  try {
    // Variables used by the function
    seidr_aggregate_param_t param;

    po::options_description umbrella("Aggregate multiple SeidrFiles.");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.out_file)->default_value("aggregated.sf"),
      "Output file name")(
      "tempdir,T",
      po::value<std::string>(&param.tempdir)->default_value("", "auto"),
      "Directory to store temporary data");

    po::options_description ompopt("OpenMP Options");
    ompopt.add_options()(
      "threads,O",
      po::value<int>(&param.nthreads)->default_value(GET_MAX_PSTL_THREADS()),
      "Number of OpenMP threads for parallel sorting");

    po::options_description agropt("Aggregate Options");
    agropt.add_options()(
      "method,m",
      po::value<std::string>(&param.method)->default_value("irp"),
      "Method to aggregate networks [top1, top2, borda, irp]")(
      "keep,k",
      po::bool_switch(&param.keep_di)->default_value(false),
      "Keep directionality information");

    po::options_description req("Required Options");
    req.add_options()("in-file",
                      po::value<std::vector<std::string>>(&param.infs),
                      "Input files");

    umbrella.add(req).add(agropt).add(ompopt).add(opt);

    po::positional_options_description p;
    p.add("in-file", -1);

    po::variables_map vm;
    po::store(
      po::command_line_parser(args).options(umbrella).positional(p).run(), vm);

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);

    set_pstl_threads(param.nthreads);
    log(LOG_INFO) << param.out_file << '\n';
    param.out_file = to_absolute(param.out_file);

    if (!param.force) {
      assert_no_overwrite(param.out_file);
    }

    assert_dir_is_writeable(dirname(param.out_file));

    for (auto& inf : param.infs) {
      inf = to_canonical(inf);
      assert_can_read(inf);
      log(LOG_INFO) << "Have file " << inf << '\n';
    }
    assert_arg_constraint<std::string>({ "top1", "top2", "borda", "irp" },
                                       param.method);

    std::vector<SeidrFile> infile_vec;
    std::vector<SeidrFileHeader> header_vec;

    log(LOG_INFO) << "Getting headers\n";
    uint32_t counter = 0;
    for (auto& f : param.infs) {
      // Create a SeidrFile from input string and call bgzf_open()
      infile_vec.emplace_back(SeidrFile(f.c_str()));
      infile_vec[counter].open("r");
      // Read file header and push into vector
      SeidrFileHeader h;
      h.unserialize(infile_vec[counter]);
      header_vec.push_back(h);
      counter++;
    }

    log(LOG_DEBUG) << "Create new header\n";
    // Setup a new header (temporary) for serializing
    SeidrFileHeader h;
    h.attr.nodes = header_vec[0].attr.nodes;
    h.attr.dense = 1;
    h.attr.nalgs = infile_vec.size();
    for (auto& hh : header_vec) {
      h.attr.nsupp += hh.attr.nsupp;
      h.attr.nsupp_str += hh.attr.nsupp_str;
      h.attr.nsupp_int += hh.attr.nsupp_int;
      h.attr.nsupp_flt += hh.attr.nsupp_flt;
      for (const auto& s : hh.supp) {
        h.supp.push_back(s);
      }
      for (const auto& a : hh.algs) {
        h.algs.push_back(a);
      }
    }
    h.nodes = header_vec[0].nodes;
    h.version_from_char(_XSTR(VERSION));
    h.cmd_from_args(args, "seidr aggregate");
    h.algs.push_back(param.method);
    h.attr.nalgs += 1;

    if (param.method == "top1") {
      h.attr.nsupp += 1;
      h.attr.nsupp_int += 1;
      h.supp.emplace_back("T");
    }

    if (param.method == "top2") {
      h.attr.nsupp += 2;
      h.attr.nsupp_int += 2;
      h.supp.emplace_back("T1");
      h.supp.emplace_back("T2");
    }

    if (param.keep_di) {
      h.attr.nsupp += h.attr.nalgs - 1;
      h.attr.nsupp_int += h.attr.nalgs - 1;
      for (uint16_t i = 1; i < h.attr.nalgs; i++) {
        h.supp.push_back("D" + std::to_string(i));
      }
    }

    log(LOG_INFO) << "Aggregating using method: " << param.method << '\n';

    param.tempfile = tempfile(param.tempdir);
    log(LOG_INFO) << "Using temp file: " << param.tempfile << '\n';

    SeidrFile tmp(param.tempfile.c_str());
    tmp.open("w");
    h.serialize(tmp);
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    // Aggregate in lower triangular order (guaranteed from `seidr import`)
    // First, setup vectors of length equal to the number of input files
    // One holds the edges in the input files, the other decides if we need
    // to read the next edge from this file
    std::vector<SeidrFileEdge> v;
    std::vector<uint8_t> get_next;
    std::vector<uint64_t> edges_read;
    v.resize(counter);
    get_next.resize(counter);
    edges_read.resize(counter);
    for (auto& e : edges_read) {
      e = 0;
    }
    // Define aggregation method
    std::function<void(std::vector<SeidrFileHeader>&,
                       std::vector<SeidrFileEdge>&,
                       SeidrFileEdge&,
                       std::vector<uint8_t>&,
                       std::vector<uint8_t>&,
                       double&,
                       uint64_t&,
                       bool&)>
      aggr_fun;
    if (param.method == "top1") {
      aggr_fun = aggr_top1;
    } else if (param.method == "borda") {
      aggr_fun = aggr_borda;
    } else if (param.method == "top2") {
      aggr_fun = aggr_top2;
    } else if (param.method == "irp") {
      aggr_fun = aggr_irp;
    } else {
      throw std::runtime_error("Unknown ranking method: " + param.method);
    }
    // In the first run, we read once from all files
    for (auto& x : get_next) {
      x = 1;
    }
    uint64_t index = 0;
    std::vector<aggr_rank_t> rvec;
    for (uint64_t i = 1; i < h.attr.nodes; i++) {
      for (uint64_t j = 0; j < i; j++) {
        for (uint32_t k = 0; k < counter; k++) {
          if (get_next[k] != 0) {
            edges_read[k]++;
            // If the last edge in a sparse network was read, do not attempt a
            // file read
            if (header_vec[k].attr.dense == 0 &&
                edges_read[k] > header_vec[k].attr.edges) {
              get_next[k] = 0;
              SeidrFileEdge e;
              e.index.i = std::numeric_limits<uint32_t>::infinity();
              e.index.j = std::numeric_limits<uint32_t>::infinity();
              v[k] = e;
            } else {
              SeidrFileEdge e;
              e.unserialize(infile_vec[k], header_vec[k]);
              if (header_vec[k].attr.dense != 0) {
                e.index.i = i;
                e.index.j = j;
              }
              v[k] = e;
            }
          }
        }
        SeidrFileEdge res =
          calc_score(header_vec, v, get_next, i, j, aggr_fun, param.keep_di);
        if (EDGE_EXISTS(res.attr.flag)) {
          h.attr.edges++;
          aggr_rank_t r;
          r.index = index;
          r.score = res.scores[counter].r;
          rvec.push_back(r);
          if (res.scores[counter].r < min) {
            min = res.scores[counter].r;
          }
          if (res.scores[counter].r > max) {
            max = res.scores[counter].r;
          }
        }
        res.serialize(tmp, h);
        index++;
        if (index % 100000 == 0) {
          log(LOG_INFO) << "Processed " << index << " edges\n";
        }
      }
    }
    log(LOG_INFO)
      << "Done aggregating. Ranking and standardizing aggregated edges\n";
    // call bgzf_close() on all open files
    for (auto& f : infile_vec) {
      f.close();
    }
    tmp.close();
    tmp.open("r");

    rank_vector(rvec);

    SeidrFileHeader ha;
    ha.unserialize(tmp);

    SeidrFile out(param.out_file.c_str());
    out.open("w");

    // Determine storage model of final file
    double ne = h.attr.edges;
    double nn = h.attr.nodes;
    if (ne < ((nn * (nn - 1) / 2) * 0.66)) {
      h.attr.dense = 0;
    }
    // Write header
    h.serialize(out);

    // Serialize all edges
    index = 0;
    for (uint64_t i = 1; i < h.attr.nodes; i++) {
      for (uint64_t j = 0; j < i; j++) {
        SeidrFileEdge e;
        e.unserialize(tmp, ha);
        if (EDGE_EXISTS(e.attr.flag)) {
          log(LOG_DEBUG) << rvec[index].index << ", " << e.scores[counter].r
                         << ", " << rvec[index].score << ", " << index << '\n';
          e.scores[counter].s = unity_stand(min, max, e.scores[counter].r);
          e.scores[counter].r = rvec[index].rank;
          index++;
        }
        if (h.attr.dense != 0) {
          e.serialize(out, h);
        } else if (EDGE_EXISTS(e.attr.flag)) {
          e.index.i = i;
          e.index.j = j;
          e.serialize(out, h);
        }
      }
    }

    tmp.close();
    out.close();

    log(LOG_INFO) << "Removing temp file: " << param.tempfile << '\n';
    fs::remove(fs::path(param.tempfile));
    log(LOG_INFO) << "Finished\n";

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

  return 0;
}
