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
#include <common.h>
#include <fs.h>
#include <node2vec.h>
#include <stats.h>

#include <BSlogger.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace po = boost::program_options;

#undef DEBUG

#include <networkit/embedding/Node2Vec.hpp>
#include <networkit/graph/Graph.hpp>

int node2vec(const std::vector<std::string>& args) {
  logger log(std::cerr, "node2vec");

  try {
    seidr_node2vec_param_t param;

    po::options_description umbrella(
        "Calculate embedding vectors for nodes in a seidr graph");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
        "help,h",
        "Show this help message")("outfile,o",
                                  po::value<std::string>(&param.outfile)
                                      ->default_value("-"),
                                  "Output file name ['-' for stdout]");

    po::options_description nopt("node2vec Options");
    nopt.add_options()("index,i",
      po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
      "Index of scores that should be used as weights")
    ("directed,d",
                      po::bool_switch(&param.directed)->default_value(false),
                      "Build directed graph")
    ("walk-return,P",
      po::value<double>(&param.walk_return)->default_value(SEIDR_NODE2VEC_DEF_WALK_RETURN, to_rounded_str(SEIDR_NODE2VEC_DEF_WALK_RETURN)),
      "Walk return parameter (stay local)")
    ("walk-inout,Q",
      po::value<double>(&param.walk_inout)->default_value(SEIDR_NODE2VEC_DEF_INOUT, to_rounded_str(SEIDR_NODE2VEC_DEF_WALK_RETURN)),
      "Walk in/out parameter (drift away)")
    ("walk-length,L",
      po::value<uint64_t>(&param.walk_length)->default_value(SEIDR_NODE2VEC_DEF_WALK_LENGTH),
      "Walk length")
    ("walk-count,N",
      po::value<uint64_t>(&param.walk_count)->default_value(SEIDR_NODE2VEC_DEF_WALK_COUNT),
      "Walk count")
    ("dimension,D",
      po::value<uint32_t>(&param.dimension)->default_value(SEIDR_NODE2VEC_DEF_DIMENSION),
      "Dimension of embedding vectors");

    po::options_description req("Required Options [can be positional]");
    req.add_options()("in-file",
                      po::value<std::string>(&param.infile)->required(),
                      "Input SeidrFile");

    umbrella.add(req).add(nopt).add(opt);

    po::positional_options_description p;
    p.add("in-file", 1);

    po::variables_map vm;
    po::store(
        po::command_line_parser(args).options(umbrella).positional(p).run(),
        vm);

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);
    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);
    if (param.outfile != "-") {
      param.outfile = to_absolute(param.outfile);
      if (!param.force) {
        assert_no_overwrite(param.outfile);
      }
      assert_dir_is_writeable(dirname(param.outfile));
    }

    std::shared_ptr<std::ostream> out;
    if (param.outfile == "-") {
      out.reset(&std::cout, no_delete);
    } else {
      out = std::shared_ptr<std::ostream>(
          new std::ofstream(param.outfile.c_str()));
    }

    SeidrFile rf(param.infile.c_str());
    rf.open("r");

    SeidrFileHeader h;
    h.unserialize(rf);

    make_tpos(param.tpos, h);

    log(LOG_INFO) << "Reading graph\n";
    std::vector<SeidrFileEdge> ev = read_network(
        h, rf, -std::numeric_limits<double>::infinity(), param.tpos, false);

    log(LOG_INFO) << "Creating Networkit graph\n";
    NetworKit::Graph g(h.attr.nodes, true, param.directed);
    for (auto& e : ev) {
      if (!param.directed) {
        g.addEdge(e.index.i, e.index.j, e.scores[param.tpos].s);
      } else {
        if (EDGE_IS_DIRECT(e.attr.flag)) {
          if (EDGE_IS_AB(e.attr.flag)) {
            g.addEdge(e.index.i, e.index.j, e.scores[param.tpos].s);
          } else {
            g.addEdge(e.index.j, e.index.i, e.scores[param.tpos].s);
          }
        } else {
          g.addEdge(e.index.i, e.index.j, e.scores[param.tpos].s);
          g.addEdge(e.index.j, e.index.i, e.scores[param.tpos].s);
        }
      }
    }

    NetworKit::Node2Vec n2v(g, param.walk_return, param.walk_inout,
                            param.walk_length, param.walk_count,
                            param.dimension);
    n2v.run();

    std::vector<std::vector<float>> const& features = n2v.getFeatures();

    for (uint64_t i = 0; i < features.size(); i++) {
      (*out) << h.nodes[i] << "\t";
      for (uint32_t j = 0; j < param.dimension; j++) {
        (*out) << features[i][j] << (j == (param.dimension - 1) ? '\n' : '\t');
      }
    }

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
