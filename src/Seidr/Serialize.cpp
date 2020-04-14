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

#include <Serialize.h>
#include <algorithm>
#include <armadillo>
#include <bgzf.h>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef MAKE_EXTERN
extern "C"
{
#endif

  void SeidrFile::open(const char* mode)
  {
    std::string m(mode);
    if (m == "w")
      m += "5";
    bgzfile = bgzf_open(filepath, m.c_str());
    _opened = 1;
  }

  void SeidrFile::close()
  {
    if (_opened && !_closed) {
      bgzf_close(bgzfile);
      _closed = 1;
      _opened = 0;
    }
  }

  void SeidrFile::seek(int64_t offset)
  {
    int64_t r = bgzf_seek(bgzfile, offset, SEEK_SET);
    if (r == -1)
      throw std::runtime_error("Error seeking in SeidrFile\n");
  }

  int64_t SeidrFile::tell() { return bgzf_tell(bgzfile); }

  void SeidrFileHeader::serialize(SeidrFile& file) const
  {
    ssize_t n = 0;
    n = bgzf_write(file.bgzfile, &attr, sizeof(header_attr));
    if (n == 0)
      throw std::runtime_error("Couldn't write SeidrFileHeader attributes\n");
    for (uint64_t i = 0; i < attr.nodes; i++) {
      header_node hn;
      const char* cn = nodes[i].c_str();
      strncpy(hn.data, cn, HEADER_NODE_SIZE);
      n = bgzf_write(file.bgzfile, &hn, sizeof(header_node));
      if (n == 0)
        throw std::runtime_error("Couldn't write SeidrFileHeader node\n");
    }
    for (uint16_t i = 0; i < attr.nalgs; i++) {
      header_alg ha;
      const char* cn = algs[i].c_str();
      strncpy(ha.data, cn, HEADER_ALG_SIZE);
      n = bgzf_write(file.bgzfile, &ha, sizeof(header_alg));
      if (n == 0)
        throw std::runtime_error("Couldn't write SeidrFileHeader "
                                 "Algorithms\n");
    }
    for (uint32_t i = 0; i < attr.nsupp; i++) {
      header_alg hs;
      const char* cn = supp[i].c_str();
      strncpy(hs.data, cn, HEADER_SUPP_SIZE);
      n = bgzf_write(file.bgzfile, &hs, sizeof(header_supp));
      if (n == 0)
        throw std::runtime_error("Couldn't write SeidrFileHeader "
                                 "supplementary information\n");
    }
    if (attr.pagerank_calc) {
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        c.data = pagerank[i];
        n = bgzf_write(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't write SeidrFileHeader "
                                   "pagerank\n");
      }
    }
    if (attr.closeness_calc) {
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        c.data = closeness[i];
        n = bgzf_write(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't write SeidrFileHeader "
                                   "closeness\n");
      }
    }
    if (attr.betweenness_calc) {
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        c.data = betweenness[i];
        n = bgzf_write(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't write SeidrFileHeader "
                                   "betweenness\n");
      }
    }
    if (attr.strength_calc) {
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        c.data = strength[i];
        n = bgzf_write(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't write SeidrFileHeader "
                                   "strength\n");
      }
    }
    if (attr.eigenvector_calc) {
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        c.data = eigenvector[i];
        n = bgzf_write(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't write SeidrFileHeader "
                                   "eigenvector\n");
      }
    }
    if (attr.katz_calc) {
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        c.data = katz[i];
        n = bgzf_write(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't write SeidrFileHeader "
                                   "katz\n");
      }
    }
    int success = bgzf_flush(file.bgzfile);
    if (success != 0)
      throw std::runtime_error("Couldn't flush SeidrFileHeader\n");
  }

  void SeidrFileHeader::unserialize(SeidrFile& file)
  {
    ssize_t n = 0;
    n = bgzf_read(file.bgzfile, &attr, sizeof(header_attr));
    if (n == 0)
      throw std::runtime_error("Couldn't read SeidrFileHeader attributes\n");
    nodes.resize(attr.nodes);
    algs.resize(attr.nalgs);
    supp.resize(attr.nsupp);
    for (uint64_t i = 0; i < attr.nodes; i++) {
      header_node hn;
      n = bgzf_read(file.bgzfile, &hn, sizeof(header_node));
      if (n == 0)
        throw std::runtime_error("Couldn't read SeidrFileHeader nodes\n");
      nodes[i] = std::string(hn.data);
    }
    for (uint64_t i = 0; i < attr.nalgs; i++) {
      header_alg ha;
      n = bgzf_read(file.bgzfile, &ha, sizeof(header_alg));
      if (n == 0)
        throw std::runtime_error("Couldn't read SeidrFileHeader algorithms\n");
      algs[i] = std::string(ha.data);
    }
    for (uint32_t i = 0; i < attr.nsupp; i++) {
      header_supp hs;
      n = bgzf_read(file.bgzfile, &hs, sizeof(header_supp));
      if (n == 0)
        throw std::runtime_error("Couldn't read SeidrFileHeader "
                                 "supplementary information\n");
      supp[i] = std::string(hs.data);
    }
    if (attr.pagerank_calc) {
      pagerank.resize(attr.nodes);
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        n = bgzf_read(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't read SeidrFileHeader "
                                   "pagerank\n");
        pagerank[i] = c.data;
      }
    }
    if (attr.closeness_calc) {
      closeness.resize(attr.nodes);
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        n = bgzf_read(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't read SeidrFileHeader "
                                   "closeness\n");
        closeness[i] = c.data;
      }
    }
    if (attr.betweenness_calc) {
      betweenness.resize(attr.nodes);
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        n = bgzf_read(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't read SeidrFileHeader "
                                   "betweenness\n");
        betweenness[i] = c.data;
      }
    }
    if (attr.strength_calc) {
      strength.resize(attr.nodes);
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        n = bgzf_read(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't read SeidrFileHeader "
                                   "strength\n");
        strength[i] = c.data;
      }
    }
    if (attr.eigenvector_calc) {
      eigenvector.resize(attr.nodes);
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        n = bgzf_read(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't read SeidrFileHeader "
                                   "eigenvector\n");
        eigenvector[i] = c.data;
      }
    }
    if (attr.katz_calc) {
      katz.resize(attr.nodes);
      for (uint32_t i = 0; i < attr.nodes; i++) {
        header_centrality c;
        n = bgzf_read(file.bgzfile, &c, sizeof(header_centrality));
        if (n == 0)
          throw std::runtime_error("Couldn't read SeidrFileHeader "
                                   "katz\n");
        katz[i] = c.data;
      }
    }
  }

  void SeidrFileHeader::print(std::ostream& out)
  {
    out << "# [G] Nodes: " << attr.nodes << '\n';
    out << "# [G] Edges: " << attr.edges << '\n';
    out << "# [G] Storage: " << (attr.dense ? "Dense\n" : "Sparse\n");
    out << "# [G] Algorithms #: " << attr.nalgs << '\n';
    out << "# [G] Supplementary data #: " << attr.nsupp << '\n';
    for (auto& a : algs)
      out << "# [A] " << a << '\n';
    for (auto& s : supp)
      out << "# [S] " << s << '\n';
    out << "# [R] Version: " << attr.version << '\n';
    out << "# [R] Cmd: " << attr.cmd << '\n';
    for (auto& n : nodes)
      out << "# [N] " << n << '\n';
  }

  void SeidrFileHeader::print_centrality(std::ostream& out) const
  {
    uint32_t nmeasures = attr.pagerank_calc + attr.closeness_calc +
                         attr.betweenness_calc + attr.strength_calc +
                         attr.eigenvector_calc + attr.katz_calc;
    if (nmeasures == 0)
      throw std::runtime_error(
        "At least one node centrality measure must be calculated\n");
    for (uint32_t i = 0; i < attr.nodes; i++) {
      auto nm = nmeasures;
      (out) << nodes[i] << '\t';
      if (attr.pagerank_calc) {
        (out) << pagerank[i];
        nm--;
      }
      if (nm == 0)
        (out) << '\n';
      else
        (out) << '\t';
      if (attr.closeness_calc) {
        (out) << closeness[i];
        nm--;
      }
      if (nm == 0)
        (out) << '\n';
      else
        (out) << '\t';
      if (attr.betweenness_calc) {
        (out) << betweenness[i];
        nm--;
      }
      if (nm == 0)
        (out) << '\n';
      else
        (out) << '\t';
      if (attr.strength_calc) {
        (out) << strength[i];
        nm--;
      }
      if (nm == 0)
        (out) << '\n';
      else
        (out) << '\t';
      if (attr.eigenvector_calc) {
        (out) << eigenvector[i];
        nm--;
      }
      if (nm == 0)
        (out) << '\n';
      else
        (out) << '\t';
      if (attr.katz_calc) {
        (out) << katz[i];
        nm--;
      }
      if (nm == 0)
        (out) << '\n';
      else
        (out) << '\t';
    }
  }

  void SeidrFileHeader::version_from_char(const char* version)
  {
    strncpy(attr.version, version, HEADER_VERSION_SIZE);
  }

  void SeidrFileHeader::cmd_from_args(char** argv, int argc)
  {
    std::string cmd = "";
    for (int i = 0; i <= argc; i++) {
      std::string tmp = argv[i];
      cmd += tmp;
      if (i < argc)
        cmd += " ";
    }
    const char* c = cmd.c_str();
    strncpy(attr.cmd, c, HEADER_CMD_SIZE);
  }

  uint16_t SeidrFileHeader::get_supp_ind(std::string supp_n)
  {
    if (attr.nsupp == 0)
      throw std::runtime_error("No supplementary data.");

    auto index_it = std::find(supp.begin(), supp.end(), supp_n);

    if (index_it == supp.end())
      throw std::runtime_error("Supplementary data: " + supp_n + " not found.");

    uint16_t offset = std::distance(supp.begin(), index_it);

    if (offset < attr.nsupp_str)
      return offset;
    else if (offset < (attr.nsupp_int + attr.nsupp_str))
      return offset - attr.nsupp_str;
    else
      return offset - attr.nsupp_str - attr.nsupp_int;
  }

  bool SeidrFileHeader::have_supp(std::string supp_n)
  {
    if (attr.nsupp == 0)
      return false;

    auto index_it = std::find(supp.begin(), supp.end(), supp_n);

    if (index_it == supp.end())
      return false;
    else
      return true;
  }

  void SeidrFileEdge::serialize(SeidrFile& file, SeidrFileHeader& header) const
  {
    ssize_t n = 0;
    if (header.attr.dense) {
      n = bgzf_write(file.bgzfile, &attr, sizeof(edge_attr));
      if (n == 0)
        raise_err("Couldn't write SeidrFileEdge attributes\n", file, header);
      if (EDGE_EXISTS(attr.flag)) {
        for (uint16_t i = 0; i < header.attr.nalgs; i++) {
          n = bgzf_write(file.bgzfile, &scores[i], sizeof(edge_score));
          if (n == 0)
            raise_err("Couldn't write SeidrFileEdge"
                      " scores\n",
                      file,
                      header);
        }
        for (uint16_t i = 0; i < header.attr.nsupp_str; i++) {
          const char* sinfo = supp_str[i].c_str();
          edge_supp_str es;
          strncpy(es.data, sinfo, EDGE_SUPP_SIZE);
          n = bgzf_write(file.bgzfile, &es, sizeof(edge_supp_str));
          if (n == 0)
            raise_err("Couldn't write SeidrFileEdge"
                      " supplementary data\n",
                      file,
                      header);
        }
        for (uint16_t i = 0; i < header.attr.nsupp_int; i++) {
          edge_supp_int es;
          es.data = supp_int[i];
          n = bgzf_write(file.bgzfile, &es, sizeof(edge_supp_int));
          if (n == 0)
            raise_err("Couldn't write SeidrFileEdge"
                      " supplementary data\n",
                      file,
                      header);
        }
        for (uint16_t i = 0; i < header.attr.nsupp_flt; i++) {
          edge_supp_flt es;
          es.data = supp_flt[i];
          n = bgzf_write(file.bgzfile, &es, sizeof(edge_supp_flt));
          if (n == 0)
            raise_err("Couldn't write SeidrFileEdge"
                      " supplementary data\n",
                      file,
                      header);
        }
      }
    } else {
      n = bgzf_write(file.bgzfile, &index, sizeof(edge_index));
      if (n == 0)
        raise_err("Couldn't write SeidrFileEdge index\n", file, header);
      n = bgzf_write(file.bgzfile, &attr, sizeof(edge_attr));
      if (n == 0)
        raise_err("Couldn't write SeidrFileEdge attributes\n", file, header);
      for (uint16_t i = 0; i < header.attr.nalgs; i++) {
        n = bgzf_write(file.bgzfile, &scores[i], sizeof(edge_score));
        if (n == 0)
          raise_err("Couldn't write SeidrFileEdge"
                    " scores\n",
                    file,
                    header);
      }
      for (uint16_t i = 0; i < header.attr.nsupp_str; i++) {
        const char* sinfo = supp_str[i].c_str();
        edge_supp_str es;
        strncpy(es.data, sinfo, EDGE_SUPP_SIZE);
        n = bgzf_write(file.bgzfile, &es, sizeof(edge_supp_str));
        if (n == 0)
          raise_err("Couldn't write SeidrFileEdge"
                    " supplementary data\n",
                    file,
                    header);
      }
      for (uint16_t i = 0; i < header.attr.nsupp_int; i++) {
        edge_supp_int es;
        es.data = supp_int[i];
        n = bgzf_write(file.bgzfile, &es, sizeof(edge_supp_int));
        if (n == 0)
          raise_err("Couldn't write SeidrFileEdge"
                    " supplementary data\n",
                    file,
                    header);
      }
      for (uint16_t i = 0; i < header.attr.nsupp_flt; i++) {
        edge_supp_flt es;
        es.data = supp_flt[i];
        n = bgzf_write(file.bgzfile, &es, sizeof(edge_supp_flt));
        if (n == 0)
          raise_err("Couldn't write SeidrFileEdge"
                    " supplementary data\n",
                    file,
                    header);
      }
    }
  }

  void SeidrFileEdge::unserialize(SeidrFile& file, SeidrFileHeader& header)
  {
    ssize_t n = 0;

    if (header.attr.dense) {
      n = bgzf_read(file.bgzfile, &attr, sizeof(edge_attr));
      if (n == 0)
        raise_err("Couldn't read SeidrFileEdge attributes\n", file, header);
      if (EDGE_EXISTS(attr.flag)) {
        scores.resize(header.attr.nalgs);
        supp_str.resize(header.attr.nsupp_str);
        supp_int.resize(header.attr.nsupp_int);
        supp_flt.resize(header.attr.nsupp_flt);
        for (uint16_t i = 0; i < header.attr.nalgs; i++) {
          n = bgzf_read(file.bgzfile, &scores[i], sizeof(edge_score));
          if (n == 0)
            raise_err("Couldn't read SeidrFileEdge"
                      " scores\n",
                      file,
                      header);
        }
        for (uint16_t i = 0; i < header.attr.nsupp_str; i++) {
          edge_supp_str es;
          n = bgzf_read(file.bgzfile, &es, sizeof(edge_supp_str));
          if (n == 0)
            raise_err("Couldn't read SeidrFileEdge"
                      " supplementary data\n",
                      file,
                      header);
          std::string sinfo(es.data);
          supp_str[i] = sinfo;
        }
        for (uint16_t i = 0; i < header.attr.nsupp_int; i++) {
          edge_supp_int es;
          n = bgzf_read(file.bgzfile, &es, sizeof(edge_supp_int));
          if (n == 0)
            raise_err("Couldn't read SeidrFileEdge"
                      " supplementary data\n",
                      file,
                      header);
          supp_int[i] = es.data;
        }
        for (uint16_t i = 0; i < header.attr.nsupp_flt; i++) {
          edge_supp_flt es;
          n = bgzf_read(file.bgzfile, &es, sizeof(edge_supp_flt));
          if (n == 0)
            raise_err("Couldn't read SeidrFileEdge"
                      " supplementary data\n",
                      file,
                      header);
          supp_flt[i] = es.data;
        }
      }
    } else {
      n = bgzf_read(file.bgzfile, &index, sizeof(edge_index));
      if (n == 0)
        raise_err("Couldn't read SeidrFileEdge index\n", file, header);
      n = bgzf_read(file.bgzfile, &attr, sizeof(edge_attr));
      if (n == 0)
        raise_err("Couldn't read SeidrFileEdge attributes\n", file, header);
      scores.resize(header.attr.nalgs);
      supp_str.resize(header.attr.nsupp_str);
      supp_int.resize(header.attr.nsupp_int);
      supp_flt.resize(header.attr.nsupp_flt);
      for (uint16_t i = 0; i < header.attr.nalgs; i++) {
        n = bgzf_read(file.bgzfile, &scores[i], sizeof(edge_score));
        if (n == 0)
          raise_err("Couldn't read SeidrFileEdge"
                    " scores\n",
                    file,
                    header);
      }
      for (uint16_t i = 0; i < header.attr.nsupp_str; i++) {
        edge_supp_str es;
        n = bgzf_read(file.bgzfile, &es, sizeof(edge_supp_str));
        if (n == 0)
          raise_err("Couldn't read SeidrFileEdge"
                    " supplementary data\n",
                    file,
                    header);
        std::string sinfo(es.data);
        supp_str[i] = sinfo;
      }
      for (uint16_t i = 0; i < header.attr.nsupp_int; i++) {
        edge_supp_int es;
        n = bgzf_read(file.bgzfile, &es, sizeof(edge_supp_int));
        if (n == 0)
          raise_err("Couldn't read SeidrFileEdge"
                    " supplementary data\n",
                    file,
                    header);
        supp_int[i] = es.data;
      }
      for (uint16_t i = 0; i < header.attr.nsupp_flt; i++) {
        edge_supp_flt es;
        n = bgzf_read(file.bgzfile, &es, sizeof(edge_supp_flt));
        if (n == 0)
          raise_err("Couldn't read SeidrFileEdge"
                    " supplementary data\n",
                    file,
                    header);
        supp_flt[i] = es.data;
      }
    }
  }

  void SeidrFileEdge::print(std::ostream& out,
                            const SeidrFileHeader& header,
                            char end,
                            bool print_supp,
                            std::string delim,
                            std::string sc_delim,
                            bool full,
                            bool supp_tag,
                            bool no_name) const
  {
    uint32_t i = index.i;
    uint32_t j = index.j;
    uint32_t xi, xj;

    if (i < j) {
      uint32_t itmp = i;
      i = j;
      j = itmp;
    }

    if (EDGE_IS_DIRECT(attr.flag)) {
      if (EDGE_IS_BA(attr.flag)) {
        if (no_name) {
          out << j << '\t' << i << '\t' << "Directed\t";
        } else {
          out << header.nodes[j] << '\t' << header.nodes[i] << '\t'
              << "Directed\t";
        }
        xi = j;
        xj = i;
      } else {
        if (no_name) {
          out << i << '\t' << j << '\t' << "Directed\t";
        } else {
          out << header.nodes[i] << '\t' << header.nodes[j] << '\t'
              << "Directed\t";
        }
        xi = i;
        xj = j;
      }
    } else {
      if (no_name) {
        out << i << '\t' << j << '\t' << "Undirected\t";
      } else {
        out << header.nodes[i] << '\t' << header.nodes[j] << '\t'
            << "Undirected\t";
      }
      xi = i;
      xj = j;
    }
    for (uint16_t k = 0; k < header.attr.nalgs; k++) {
      out << scores[k].s << delim << scores[k].r;
      if (k != header.attr.nalgs - 1)
        out << '\t';
    }
    if (print_supp && header.attr.nsupp > 0) {
      out << '\t';
      uint32_t head_inx = 0;
      for (uint16_t k = 0; k < header.attr.nsupp_str; k++) {
        if (supp_tag)
          out << header.supp[head_inx] << ":";
        out << supp_str[k];
        if (head_inx < header.attr.nsupp - 1)
          out << sc_delim;
        head_inx++;
      }
      for (uint16_t k = 0; k < header.attr.nsupp_int; k++) {
        if (supp_tag)
          out << header.supp[head_inx] << ":";
        out << supp_int[k];
        if (head_inx < header.attr.nsupp - 1)
          out << sc_delim;
        head_inx++;
      }
      for (uint16_t k = 0; k < header.attr.nsupp_flt; k++) {
        if (supp_tag)
          out << header.supp[head_inx] << ":";
        out << supp_flt[k];
        if (head_inx < header.attr.nsupp - 1)
          out << sc_delim;
        head_inx++;
      }
    }
    if (!full) {
      out << end;
    } else {
      std::stringstream ost;
      ost << '\t';
      if (header.attr.pagerank_calc)
        ost << header.pagerank[xi] << '\t' << header.pagerank[xj] << '\t';
      if (header.attr.closeness_calc)
        ost << header.closeness[xi] << '\t' << header.closeness[xj] << '\t';
      if (header.attr.betweenness_calc)
        ost << header.betweenness[xi] << '\t' << header.betweenness[xj] << '\t';
      if (header.attr.strength_calc)
        ost << header.strength[xi] << '\t' << header.strength[xj] << '\t';
      if (header.attr.eigenvector_calc)
        ost << header.eigenvector[xi] << '\t' << header.eigenvector[xj] << '\t';
      if (header.attr.katz_calc)
        ost << header.katz[xi] << '\t' << header.katz[xj] << '\t';
      std::string ox = ost.str();
      if (ox.size() > 0) {
        ox.back() = '\n';
      } else {
        ox += '\n';
      }
      out << ox;
    }
  }

  void SeidrFileEdge::raise_err(std::string what,
                                SeidrFile& file,
                                SeidrFileHeader& header) const
  {
    std::string f(file.filepath);
    std::string err = "File '" + f + "', at [" + std::to_string(index.i) + "," +
                      std::to_string(index.j) + "]: " + what;
    throw std::runtime_error(err);
  }

  void SeidrFileIndex::build(SeidrFile& inf)
  {
    SeidrFileHeader h;
    h.unserialize(inf);
    data.resize((h.attr.nodes * (h.attr.nodes - 1)) / 2, 0);
    nodes = h.attr.nodes;
    edges = h.attr.edges;
    if (h.attr.dense) {
      dense = 1;
      uint64_t ctr = 0;
      for (uint32_t i = 1; i < h.attr.nodes; i++) {
        for (uint32_t j = 0; j < i; j++) {
          int64_t offset = bgzf_tell(inf.bgzfile);
          SeidrFileEdge e;
          e.unserialize(inf, h);
          if (EDGE_EXISTS(e.attr.flag)) {
            data[ctr] = offset;
          }
          ctr++;
        }
      }
    } else {
      dense = 0;
      for (uint64_t i = 0; i < h.attr.edges; i++) {
        int64_t offset = bgzf_tell(inf.bgzfile);
        SeidrFileEdge e;
        e.unserialize(inf, h);
        uint64_t inx = coord_to_inx(e.index.i, e.index.j);
        data[inx] = offset;
      }
    }
  }

  void SeidrFileIndex::build(std::string& inf)
  {
    SeidrFile f(inf.c_str());
    f.open("r");
    build(f);
    f.close();
  }

  void SeidrFileIndex::serialize(SeidrFile& file)
  {
    ssize_t n = 0;
    index_header h;
    h.size = data.size();
    h.nodes = nodes;
    h.edges = edges;
    if (dense) {
      h.dense = 1;
      n = bgzf_write(file.bgzfile, &h, sizeof(index_header));
      //
      for (uint64_t i = 0; i < data.size(); i++) {
        index_offset off;
        off.data = data[i];
        n = bgzf_write(file.bgzfile, &off, sizeof(index_offset));
      }
    } else {
      h.dense = 0;
      n = bgzf_write(file.bgzfile, &h, sizeof(index_header));
      //
      uint64_t ctr = 0;
      for (uint32_t i = 1; i < nodes; i++) {
        for (uint32_t j = 0; j < i; j++) {
          index_index coord;
          coord.i = i;
          coord.j = j;
          if (data[ctr]) {
            index_offset off;
            off.data = data[ctr];
            n = bgzf_write(file.bgzfile, &coord, sizeof(index_index));
            // std::cerr << sizeof(index_index) << ", " << n << '\n';
            n = bgzf_write(file.bgzfile, &off, sizeof(index_offset));
            // std::cerr << sizeof(index_offset) << ", " << n << '\n';
          }
          ctr++;
        }
      }
    }
  }

  void SeidrFileIndex::unserialize(SeidrFile& file, SeidrFileHeader& head)
  {
    ssize_t n = 0;
    index_header h;
    n = bgzf_read(file.bgzfile, &h, sizeof(index_header));
    dense = h.dense;
    nodes = h.nodes;
    edges = h.edges;
    uint32_t gmap_i = 0;
    for (const auto& node : head.nodes) {
      std::string xn = node;
      if (_case_insensitive) {
        std::for_each(
          xn.begin(), xn.end(), [](char& c) { c = std::tolower(c); });
      }
      gmap[xn] = gmap_i++;
    }
    uint64_t size = h.size;
    data.resize(size, -99);
    if (dense) {
      for (uint64_t i = 0; i < size; i++) {
        index_offset off;
        n = bgzf_read(file.bgzfile, &off, sizeof(index_offset));
        data[i] = off.data;
      }
    } else {
      for (uint64_t i = 0; i < edges; i++) {
        index_index coord;
        n = bgzf_read(file.bgzfile, &coord, sizeof(index_index));
        index_offset off;
        n = bgzf_read(file.bgzfile, &off, sizeof(index_offset));
        uint64_t inx = coord_to_inx(coord.i, coord.j);
        data[inx] = off.data;
      }
    }
  }

  std::pair<bool, uint32_t> SeidrFileIndex::find(const std::string& n)
  {
    uint32_t r = 0;
    std::string xn = n;
    if (_case_insensitive) {
      std::for_each(xn.begin(), xn.end(), [](char& c) { c = std::tolower(c); });
    }
    auto ptr = gmap.find(xn);
    if (ptr == gmap.end()) {
      if (_strict) {
        throw SeidrFileIndexNodeNotFound("Index not found for: " + xn);
      }
      return std::pair<bool, uint32_t>({ false, r });
    }
    r = ptr->second;
    return std::pair<bool, uint32_t>({ true, r });
  }

  std::vector<offset_t> SeidrFileIndex::get_offset_node(std::string& node)
  {
    std::vector<offset_t> ret;
    std::pair<bool, uint32_t> promise_n = find(node);
    if (promise_n.first) {
      uint32_t n = promise_n.second;
      for (uint32_t i = 0; i < nodes; i++) {
        if (i != n) {
          uint32_t xi = n;
          uint32_t xj = i;
          if (xi < xj) {
            uint32_t tmp = xi;
            xi = xj;
            xj = tmp;
          }
          uint64_t inx = coord_to_inx(xi, xj);
          if (data[inx]) {
            offset_t o;
            o.i = xi;
            o.j = xj;
            o.o = data[inx];
            ret.push_back(o);
          }
        }
      }
      std::sort(ret.begin(), ret.end());
    } else {
      offset_t o;
      o.err += "Node " + node + " not in header;";
      o.o = -1;
      ret.push_back(o);
    }
    return ret;
  }

  offset_t SeidrFileIndex::get_offset_pair(std::string& lhs, std::string& rhs)
  {
    offset_t o;
    std::pair<bool, uint32_t> promise_i = find(lhs);
    std::pair<bool, uint32_t> promise_j = find(rhs);
    if (promise_i.first && promise_j.first) {
      uint32_t i = promise_i.second;
      uint32_t j = promise_j.second;
      if (i == j) {
        o.o = -1;
        o.err += "Self links not allowed for " + lhs + ":" + rhs + ";";
      } else {
        if (i < j) {
          uint32_t tmp = i;
          i = j;
          j = tmp;
        }
        uint64_t inx = coord_to_inx(i, j);

        o.i = i;
        o.j = j;
        o.o = data[inx];
      }
    } else {
      if ((!promise_i.first) && promise_j.first) {
        o.err +=
          "Node '" + lhs + "' of pair '" + lhs + ":" + rhs + "' not found;";
      } else if ((!promise_j.first) && promise_i.first) {
        o.err +=
          "Node '" + rhs + "' of pair '" + lhs + ":" + rhs + "' not found";
      } else {
        o.err += "Both nodes of pair '" + lhs + ":" + rhs + "' not found";
      }
      o.o = -1;
    }
    return o;
  }

  std::set<offset_t> SeidrFileIndex::get_offset_nodelist(
    std::vector<std::string>& nodelist)
  {
    std::set<offset_t> ret;
    for (auto& n : nodelist) {
      bool have_pair = false;
      std::string lhs, rhs;
      for (const char& c : n) {
        if (have_pair) {
          rhs += c;
        } else if (c == ':') {
          have_pair = true;
        } else {
          lhs += c;
        }
      }
      if (have_pair) {
        offset_t off = get_offset_pair(lhs, rhs);
        ret.insert(off);
      } else {
        std::vector<offset_t> off = get_offset_node(n);
        for (auto& o : off) {
          ret.insert(o);
        }
      }
    }
    return ret;
  }

  void SeidrFile::each_edge(
    std::function<void(SeidrFileEdge&, SeidrFileHeader&)> f,
    bool include_missing)
  {
    SeidrFileHeader h;
    h.unserialize((*this));
    if (h.attr.dense) {
      for (uint64_t i = 1; i < h.attr.nodes; i++) {
        for (uint64_t j = 0; j < i; j++) {
          SeidrFileEdge e;
          e.unserialize((*this), h);
          e.index.i = i;
          e.index.j = j;
          if (EDGE_EXISTS(e.attr.flag) || include_missing) {
            f(e, h);
          }
        }
      }
    } else {
      for (uint64_t i = 0; i < h.attr.edges; i++) {
        SeidrFileEdge e;
        e.unserialize((*this), h);
        f(e, h);
      }
    }
  }

  void SeidrFile::each_edge_exit_early(
    std::function<bool(SeidrFileEdge&, SeidrFileHeader&)> f)
  {
    SeidrFileHeader h;
    h.unserialize((*this));
    if (h.attr.dense) {
      for (uint64_t i = 1; i < h.attr.nodes; i++) {
        for (uint64_t j = 0; j < i; j++) {
          SeidrFileEdge e;
          e.unserialize((*this), h);
          e.index.i = i;
          e.index.j = j;
          if (EDGE_EXISTS(e.attr.flag)) {
            if (f(e, h)) {
              return;
            }
          }
        }
      }
    } else {
      for (uint64_t i = 0; i < h.attr.edges; i++) {
        SeidrFileEdge e;
        e.unserialize((*this), h);
        if (f(e, h)) {
          return;
        }
      }
    }
  }

  std::vector<SeidrFileEdge> read_network(SeidrFileHeader& h,
                                          SeidrFile& rf,
                                          double threshold,
                                          uint32_t tpos,
                                          bool trank)
  {
    std::vector<SeidrFileEdge> ret;

    if (h.attr.dense) {
      for (uint64_t i = 1; i < h.attr.nodes; i++) {
        for (uint64_t j = 0; j < i; j++) {
          SeidrFileEdge e;
          e.unserialize(rf, h);
          e.index.i = i;
          e.index.j = j;
          if (EDGE_EXISTS(e.attr.flag)) {
            double x = trank ? e.scores[tpos].r : e.scores[tpos].s;
            bool cut = (trank ? x < threshold : x > threshold);
            if (threshold > 0)
              cut = cut && !std::isnan(x);
            if (cut)
              ret.push_back(e);
          }
        }
      }
    } else {
      for (uint64_t i = 0; i < h.attr.edges; i++) {
        SeidrFileEdge e;
        e.unserialize(rf, h);
        double x = trank ? e.scores[tpos].r : e.scores[tpos].s;
        bool cut = (trank ? x < threshold : x > threshold);
        if (threshold > 0)
          cut = cut && !std::isnan(x);
        if (cut)
          ret.push_back(e);
      }
    }
    return ret;
  }

  std::vector<MiniEdge> read_network_minimal(SeidrFileHeader& h,
                                             SeidrFile& rf,
                                             double threshold,
                                             uint32_t tpos,
                                             bool trank)
  {
    std::vector<MiniEdge> ret;

    if (h.attr.dense) {
      for (uint64_t i = 1; i < h.attr.nodes; i++) {
        for (uint64_t j = 0; j < i; j++) {
          SeidrFileEdge e;
          e.unserialize(rf, h);
          e.index.i = i;
          e.index.j = j;
          if (EDGE_EXISTS(e.attr.flag)) {
            double x = trank ? e.scores[tpos].r : e.scores[tpos].s;
            bool cut = (trank ? x < threshold : x > threshold);
            if (threshold > 0)
              cut = cut && !std::isnan(x);
            if (cut) {
              MiniEdge me;
              me.i = e.index.i;
              me.j = e.index.j;
              me.s = x;
              ret.push_back(me);
            }
          }
        }
      }
    } else {
      for (uint64_t i = 0; i < h.attr.edges; i++) {
        SeidrFileEdge e;
        e.unserialize(rf, h);
        double x = trank ? e.scores[tpos].r : e.scores[tpos].s;
        bool cut = (trank ? x < threshold : x > threshold);
        if (threshold > 0)
          cut = cut && !std::isnan(x);
        if (cut) {
          MiniEdge me;
          me.i = e.index.i;
          me.j = e.index.j;
          me.s = x;
          ret.push_back(me);
        }
      }
    }
    return ret;
  }

  arma::mat read_network_arma(SeidrFileHeader& h,
                              SeidrFile& rf,
                              double threshold,
                              uint32_t tpos,
                              bool trank)
  {
    arma::mat mat(h.nodes.size(), h.nodes.size());
    if (trank)
      mat.fill(std::numeric_limits<double>::infinity());
    else
      mat.fill(-std::numeric_limits<double>::infinity());
    if (h.attr.dense) {
      for (uint64_t i = 1; i < h.attr.nodes; i++) {
        for (uint64_t j = 0; j < i; j++) {
          SeidrFileEdge e;
          e.unserialize(rf, h);
          e.index.i = i;
          e.index.j = j;
          if (EDGE_EXISTS(e.attr.flag)) {
            double x = trank ? e.scores[tpos].r : e.scores[tpos].s;
            bool cut = (trank ? x < threshold : x > threshold);
            if (threshold > 0)
              cut = cut && !std::isnan(x);
            if (cut) {
              mat(e.index.i, e.index.j) = x;
              mat(e.index.j, e.index.i) = x;
            }
          }
        }
      }
    } else {
      for (uint64_t i = 0; i < h.attr.edges; i++) {
        SeidrFileEdge e;
        e.unserialize(rf, h);
        double x = trank ? e.scores[tpos].r : e.scores[tpos].s;
        bool cut = (trank ? x < threshold : x > threshold);
        if (threshold > 0)
          cut = cut && !std::isnan(x);
        if (cut) {
          mat(e.index.i, e.index.j) = x;
          mat(e.index.j, e.index.i) = x;
        }
      }
    }
    return mat;
  }

  bool me_score_sort(MiniEdge a, MiniEdge b) { return (a.s > b.s); }

  bool me_rank_sort(MiniEdge a, MiniEdge b) { return (a.s < b.s); }

  bool me_score_sort_abs(MiniEdge a, MiniEdge b)
  {
    return fabs(a.s) > fabs(b.s);
  }

  bool me_score_sort_rev(MiniEdge a, MiniEdge b) { return (a.s < b.s); }

  bool sfe_score_sort(SeidrFileEdge a, SeidrFileEdge b)
  {
    return (a.scores[a.scores.size() - 1].s > b.scores[b.scores.size() - 1].s);
  }

  bool sfe_rank_sort(SeidrFileEdge a, SeidrFileEdge b)
  {
    return (a.scores[a.scores.size() - 1].r < b.scores[b.scores.size() - 1].r);
  }

#ifdef MAKE_EXTERN
}
#endif
