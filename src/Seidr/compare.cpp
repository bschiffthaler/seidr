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
#include <compare.h>
#include <fs.h>
// External
#include <boost/numeric/conversion/cast.hpp>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

cindex::cindex(SeidrFileHeader& a, SeidrFileHeader& b)
{
  std::set<std::string> set_a;
  uint32_t ctr = 0;
  for (auto& gene : a.nodes) {
    map_a[gene] = ctr++;
    set_a.insert(gene);
    uni.insert(gene);
  }
  // reset counter
  ctr = 0;
  for (auto& gene : b.nodes) {
    map_b[gene] = ctr++;
    if (set_a.find(gene) != set_a.end()) {
      intersect.insert(gene);
    } else {
      diff_b.insert(gene);
    }
    uni.insert(gene);
  }
  for (auto& gene : a.nodes) {
    if (intersect.find(gene) == intersect.end()) {
      diff_a.insert(gene);
    }
  }
  ctr = 0;
  for (const auto& gene : uni) {
    map_c[gene] = ctr++;
  }
  // Create vector of unique merged node names
  for (const auto& n : uni) {
    node_uni.push_back(n);
  }
}

cindex::cindex(SeidrFileHeader& a, SeidrFileHeader& b, std::string& file_orth)
{
  orth_mode = true;
  std::ifstream ifs_orth(file_orth, std::ios::in);
  std::string line;
  while (std::getline(ifs_orth, line)) {
    std::string left;
    std::string right;
    std::stringstream ss(line);
    ss >> left;
    ss >> right;
    dict_orth[left] = right;
    dict_orth[right] = left;
  }

  std::set<std::string> set_a;
  uint32_t ctr = 0;
  for (auto& gene : a.nodes) {
    map_a[gene] = ctr++;
    set_a.insert(gene);
    uni.insert(gene);
  }
  // reset counter
  std::set<std::string> set_b;
  ctr = 0;
  for (auto& gene : b.nodes) {
    map_b[gene] = ctr++;
    set_b.insert(gene);
    uni.insert(gene);
  }

  ctr = 0;
  for (const auto& gene : uni) {
    auto ptr = dict_orth.find(gene);
    if (ptr != dict_orth.end()) {
      if (map_c.find(ptr->first) == map_c.end()) {
        map_c[ptr->first] = ctr;
        map_c[ptr->second] = ctr;
        dict_index[ctr] = string_p(ptr->first, ptr->second);
        ctr++;
      }
    } else {
      map_c[gene] = ctr;
      ctr++;
    }
  }
  // Create vector of unique merged node names
  node_uni.resize(ctr);
  for (auto& n : map_c) {
    if (dict_index.find(n.second) != dict_index.end()) {
      std::string s1 = dict_index[n.second].first;
      std::string s2 = dict_index[n.second].second;
      std::string s1s2;
      s1s2.append(s1).append(",").append(s2);
      node_uni[n.second] = s1s2;
    } else {
      node_uni[n.second] = n.first;
    }
  }
}

uint32_p
remap_index(SeidrFileEdge& e, SeidrFileHeader& h, cindex& index)
{
  // Need to copy from SeidrFileEdge& e since we can't refer to fields of a
  // packed struct.
  uint32_t i = e.index.i;
  uint32_t j = e.index.j;
  uint32_p p1(i, j);
  uint32_p p1r = index.reindex(h, e);
  return p1r;
}

SeidrFileEdge
create_edge(SeidrFileEdge& e1,
            SeidrFileHeader& h1,
            SeidrFileEdge& e2,
            SeidrFileHeader& h2,
            uint32_p& p,
            uint8_t score,
            uint32_t& ti1,
            uint32_t& ti2)
{
  SeidrFileEdge e;
  e.index.i = p.first;
  e.index.j = p.second;
  e.supp_int.push_back(score);
  e.supp_flt = { 0, 0 };
  edge_score es;
  if (score == 0) {
    es.r = e1.scores[ti1].r + e2.scores[ti2].r;
    e.supp_flt[0] = boost::numeric_cast<float>(e1.scores[ti1].r);
    e.supp_flt[1] = boost::numeric_cast<float>(e2.scores[ti2].r);
  } else if (score == 1) {
    es.r = e1.scores[ti1].r;
    e.supp_flt[0] = boost::numeric_cast<float>(e1.scores[ti1].r);
  } else {
    es.r = e2.scores[ti2].r;
    e.supp_flt[1] = boost::numeric_cast<float>(e2.scores[ti2].r);
  }
  EDGE_SET_EXISTING(e.attr.flag);
  e.scores.push_back(es);
  return e;
}

void
finalize(seidr_compare_param_t& param,
         SeidrFileHeader& h,
         double min,
         double max)
{
  SeidrFile out(param.tempfile.c_str());
  out.open("r");
  SeidrFileHeader ho;
  ho.unserialize(out);

  if (param.file_out == "-") {
    if (h.attr.dense != 0) {
      throw std::runtime_error("Serializing dense output is not implemented");
    }
    for (uint64_t i = 0; i < h.attr.edges; i++) {
      SeidrFileEdge e;
      e.unserialize(out, ho);
      e.scores[0].s = unity_stand(min, max, e.scores[0].r);
      e.print(h);
    }
  } else {
    SeidrFile fin(param.file_out.c_str());
    fin.open("w");
    h.serialize(fin);

    if (h.attr.dense != 0) {
      throw std::runtime_error("Serializing dense output is not implemented");
    }
    for (uint64_t i = 0; i < h.attr.edges; i++) {
      SeidrFileEdge e;
      e.unserialize(out, ho);
      e.scores[0].s = unity_stand(min, max, e.scores[0].r);
      e.serialize(fin, h);
    }

    fin.close();
  }
  out.close();
  remove(param.tempfile);
}

void
compare_nodes(SeidrFile& f1,
              SeidrFileHeader& h1,
              SeidrFile& f2,
              SeidrFileHeader& h2,
              seidr_compare_param_t& param)
{
  cindex index;
  if (param.trans.empty()) {
    index = cindex(h1, h2);
  } else {
    index = cindex(h1, h2, param.trans);
  }

  std::shared_ptr<std::ostream> out;
  if (param.file_out == "-") {
    out.reset(&std::cout, no_delete);
  } else {
    out =
      std::shared_ptr<std::ostream>(new std::ofstream(param.file_out.c_str()));
  }

  std::vector<std::string> node_names = index.node_union();
  std::map<std::string, int> nm;
  for (std::string& n : node_names) {
    nm[n] = 3;
  }
  f1.seek(0);
  f1.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h) {
    uint32_p re = index.reindex(h, e);
    nm[node_names[re.first]] = 1;
    nm[node_names[re.second]] = 1;
  });
  f2.seek(0);
  f2.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h) {
    uint32_p re = index.reindex(h, e);
    auto ptr = nm.find(node_names[re.first]);
    if (ptr->second == 1 || ptr->second == 0) {
      ptr->second = 0;
    } else {
      ptr->second = 2;
    }

    ptr = nm.find(node_names[re.second]);
    if (ptr->second == 1 || ptr->second == 0) {
      ptr->second = 0;
    } else {
      ptr->second = 2;
    }
  });
  for (auto& n : nm) {
    (*out) << n.first << '\t' << n.second << '\n';
  }
  f1.close();
  f2.close();
}

void
next_edge(SeidrFileEdge& e,
          SeidrFile& f,
          SeidrFileHeader& h,
          uint64_t& i,
          uint64_t& j)
{
  if (h.attr.dense != 0) {
    e.unserialize(f, h);
    e.index.i = i;
    e.index.j = j;
    if (j == (i - 1)) {
      i++;
      j = 0;
    } else {
      j++;
    }
    if (!EDGE_EXISTS(e.attr.flag)) {
      next_edge(e, f, h, i, j);
    }
  } else {
    e.unserialize(f, h);
  }
}

void
compare_ss(SeidrFile& f1,
           SeidrFileHeader& h1,
           SeidrFile& f2,
           SeidrFileHeader& h2,
           SeidrFileHeader& h3,
           seidr_compare_param_t& param)
{

  SeidrFile out(param.tempfile.c_str());
  out.open("w");

  uint64_t e1i = 1;
  uint64_t e1j = 0;
  uint64_t e2i = 1;
  uint64_t e2j = 0;

  cindex index;
  if (param.trans.empty()) {
    index = cindex(h1, h2);
  } else {
    index = cindex(h1, h2, param.trans);
  }

  SeidrFileEdge e1;
  next_edge(e1, f1, h1, e1i, e1j);

  SeidrFileEdge e2;
  next_edge(e2, f2, h2, e2i, e2j);

  uint32_p p1r = remap_index(e1, h1, index);
  uint32_p p2r = remap_index(e2, h2, index);

  uint64_t ctr1 = 1;
  uint64_t ctr2 = 1;
  uint64_t ctr3 = 0;

  bool break1 = false;
  bool break2 = false;

  std::vector<std::string>& node_names = index.node_union();

  h3.attr.nodes = node_names.size();
  h3.attr.edges = 8; // TODO(bs): Why is this 8
  h3.attr.nalgs = 1;
  h3.attr.dense = 0;
  h3.attr.nsupp = 3;
  h3.attr.nsupp_flt = 2;
  h3.attr.nsupp_int = 1;
  h3.supp = { "N", "A", "B" };
  h3.algs = { "Merge" };

  h3.nodes = node_names;

  h3.serialize(out);

  double min = std::numeric_limits<double>::infinity();
  double max = -std::numeric_limits<double>::infinity();

  while (true) {
    SeidrFileEdge e3;
    if (p1r.first < p2r.first) {
      e3 = create_edge(e1, h1, e2, h2, p1r, 1, param.tindex1, param.tindex2);
      e3.serialize(out, h3);
      // Check if this was the last edge in the file
      if (break1) {
        break;
      }
      next_edge(e1, f1, h1, e1i, e1j);
      p1r = remap_index(e1, h1, index);
      ctr1++;
      if (ctr1 == h1.attr.edges) {
        break1 = true;
      }
    } else if (p1r.first > p2r.first) {
      e3 = create_edge(e1, h1, e2, h2, p2r, 2, param.tindex1, param.tindex2);
      e3.serialize(out, h3);
      // Check if this was the last edge in the file
      if (break2) {
        break;
      }
      next_edge(e2, f2, h2, e2i, e2j);
      p2r = remap_index(e2, h2, index);
      ctr2++;
      if (ctr2 == h2.attr.edges) {
        break2 = true;
      }
    } else {
      if (p1r.second < p2r.second) {
        e3 = create_edge(e1, h1, e2, h2, p1r, 1, param.tindex1, param.tindex2);
        e3.serialize(out, h3);
        // Check if this was the last edge in the file
        if (break1) {
          break;
        }
        next_edge(e1, f1, h1, e1i, e1j);
        p1r = remap_index(e1, h1, index);
        ctr1++;
        if (ctr1 == h1.attr.edges) {
          break1 = true;
        }
      } else if (p1r.second > p2r.second) {
        e3 = create_edge(e1, h1, e2, h2, p2r, 2, param.tindex1, param.tindex2);
        e3.serialize(out, h3);
        if (break2) {
          break;
        }
        next_edge(e2, f2, h2, e2i, e2j);
        p2r = remap_index(e2, h2, index);
        ctr2++;
        if (ctr2 == h2.attr.edges) {
          break2 = true;
        }
      } else {
        e3 = create_edge(e1, h1, e2, h2, p1r, 0, param.tindex1, param.tindex2);
        e3.serialize(out, h3);
        if (break1 && break2) {
          break;
        }
        if (break1) {
          next_edge(e2, f2, h2, e2i, e2j);
          p2r = remap_index(e2, h2, index);
          ctr2++;
          ctr3++;
          break;
        }
        if (break2) {
          next_edge(e1, f1, h1, e1i, e1j);
          p1r = remap_index(e1, h1, index);
          ctr1++;
          ctr3++;
          break;
        }

        next_edge(e1, f1, h1, e1i, e1j);
        p1r = remap_index(e1, h1, index);
        next_edge(e2, f2, h2, e2i, e2j);
        p2r = remap_index(e2, h2, index);

        ctr1++;
        ctr2++;
        ctr3++;
        if (ctr1 == h1.attr.edges) {
          break1 = true;
        }
        if (ctr2 == h2.attr.edges) {
          break2 = true;
        }
      }
    }
    min = e3.scores[0].r < min ? e3.scores[0].r : min;
    max = e3.scores[0].r > max ? e3.scores[0].r : max;
  }
  SeidrFileEdge e3;
  if (break1 && !break2) {
    for (uint64_t i = ctr2 - 1; i < h2.attr.edges; i++) {
      e3 = create_edge(e1, h1, e2, h2, p2r, 2, param.tindex1, param.tindex2);
      e3.serialize(out, h3);
      if (i < h2.attr.edges - 1) {
        next_edge(e2, f2, h2, e2i, e2j);
        p2r = remap_index(e2, h2, index);
      }
      min = e3.scores[0].r < min ? e3.scores[0].r : min;
      max = e3.scores[0].r > max ? e3.scores[0].r : max;
      ctr2++;
    }
  } else if (break2 && !break1) {
    for (uint64_t i = ctr1 - 1; i < h1.attr.edges; i++) {
      e3 = create_edge(e1, h1, e2, h2, p1r, 1, param.tindex1, param.tindex2);
      e3.serialize(out, h3);
      if (i < h1.attr.edges - 1) {
        next_edge(e1, f1, h1, e1i, e1j);
        p1r = remap_index(e1, h1, index);
      }
      min = e3.scores[0].r < min ? e3.scores[0].r : min;
      max = e3.scores[0].r > max ? e3.scores[0].r : max;
      ctr1++;
    }
  }
  out.close();
  // Calc standardized scores
  uint64_t ne = ctr1 + ctr2 - ctr3 - 1;
  h3.attr.edges = ne;
  finalize(param, h3, min, max);
}

int
compare(const std::vector<std::string>& args)
{

  logger log(std::cerr, "compare");

  try {

    // Variables used by the function
    seidr_compare_param_t param;

    po::options_description umbrella("Compare edges or nodes in two networks.");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.file_out)->default_value("-"),
      "Output file name ['-' for stdout]")(
      "tempdir,T",
      po::value<std::string>(&param.tempdir)->default_value("", "auto"),
      "Directory to store temporary data");

    po::options_description copt("Compare Options");
    copt.add_options()(
      "index-a,i",
      po::value<uint32_t>(&param.tindex1)->default_value(0, "last score"),
      "Merge scores on this index for network A")(
      "index-b,j",
      po::value<uint32_t>(&param.tindex1)->default_value(0, "last score"),
      "Merge scores on this index for network B")(
      "translate,t",
      po::value<std::string>(&param.trans)->default_value(""),
      "Translate node names in network A according to this table")(
      "nodes,n",
      po::bool_switch(&param.node_mode)->default_value(false),
      "Print overlap of nodes instead of edges");

    po::options_description req("Required Options [can be positional]");
    req.add_options()("network-1",
                      po::value<std::string>(&param.net1)->required(),
                      "Input SeidrFile for network A")(
      "network-2",
      po::value<std::string>(&param.net2)->required(),
      "Input SeidrFile for network B");

    umbrella.add(req).add(copt).add(opt);

    po::positional_options_description p;
    p.add("network-1", 1);
    p.add("network-2", 1);

    po::variables_map vm;
    po::store(
      po::command_line_parser(args).options(umbrella).positional(p).run(), vm);

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return EINVAL;
    }

    po::notify(vm);

    param.net1 = to_absolute(param.net1);
    param.net2 = to_absolute(param.net2);

    param.tempfile = tempfile(param.tempdir);

    assert_exists(param.net1);
    assert_exists(param.net2);
    assert_can_read(param.net1);
    assert_can_read(param.net2);
    if (!param.force && param.file_out != "-") {
      assert_no_overwrite(param.file_out);
    }
    assert_dir_is_writeable(dirname(param.tempdir));

    SeidrFile n1(param.net1.c_str());
    n1.open("r");

    SeidrFile n2(param.net2.c_str());
    n2.open("r");

    SeidrFileHeader h1;
    h1.unserialize(n1);

    SeidrFileHeader h2;
    h2.unserialize(n2);

    make_tpos(param.tindex1, h1);

    make_tpos(param.tindex2, h2);

    if (param.node_mode) {
      compare_nodes(n1, h1, n2, h2, param);
    } else {
      SeidrFileHeader h3;
      h3.version_from_char(_XSTR(VERSION));
      h3.cmd_from_args(args, "backbone");
      compare_ss(n1, h1, n2, h2, h3, param);
    }

    n1.close();
    n2.close();

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
