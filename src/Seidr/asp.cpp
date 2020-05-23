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

#include <BSlogger.hpp>
#include <Serialize.h>
#include <asp.h>
#include <common.h>
#include <fs.h>

#undef DEBUG

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <networkit/distance/APSP.h>
#include <networkit/graph/Graph.h>
#include <string>
#include <vector>

namespace po = boost::program_options;

int
asp(const std::vector<std::string>& args)
{
  logger log(std::cerr, "asp");

  try {

    seidr_asp_param_t param;

    po::options_description umbrella(
      "Compute all shortest paths in the network");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("-"),
      "Output file name ['-' for stdout]");

    po::options_description aspopt("ASP Options");
    aspopt.add_options()(
      "index,i",
      po::value<uint32_t>(&param.tpos)->default_value(0, "last column"),
      "Algorithm position in SeidrFile to use as weights")(
      "use-rank,r",
      po::bool_switch(&param.trank)->default_value(false),
      "Use rank as edge weight basis instead of score")(
      "invert,I",
      po::bool_switch(&param.invert)->default_value(false),
      "Invert scores (if higher scores are better)")(
      "absolute,a",
      po::bool_switch(&param.absolute)->default_value(false),
      "Use absolute weights.");

    po::options_description fopt("Formatting Options");
    fopt.add_options()(
      "precision,p",
      po::value<uint16_t>(&param.precision)->default_value(ASP_DEF_PRECISION),
      "Number of decimals to print");

    po::options_description req("Required Options [can be positional]");
    req.add_options()(
      "in-file", po::value<std::string>(&param.infile), "Input SeidrFile");

    umbrella.add(req).add(aspopt).add(fopt).add(opt);

    po::positional_options_description p;
    p.add("in-file", 1);

    po::variables_map vm;
    po::store(
      po::command_line_parser(args).options(umbrella).positional(p).run(), vm);

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);

    assert_exists(param.infile);
    assert_can_read(param.infile);

    SeidrFile rf(param.infile.c_str());
    rf.open("r");

    SeidrFileHeader h;
    h.unserialize(rf);

    make_tpos(param.tpos, h);

    log(LOG_INFO) << "Starting analysis\n";

    log(LOG_INFO) << "Reading network\n";
    std::vector<MiniEdge> edges = read_network_minimal(
      h, rf, -std::numeric_limits<double>::infinity(), param.tpos, param.trank);
    log(LOG_INFO) << "Read " << edges.size() << " edges\n";

    NetworKit::Graph g(h.attr.nodes, true, false);
    for (const auto& edge : edges) {
      double weight = edge.s;
      if (almost_equal(weight, 0)) {
        weight = ASP_DEF_MIN_WEIGHT;
      }
      weight = param.absolute ? fabs(weight) : weight;
      weight = param.invert ? 1.0 / weight : weight;
      g.addEdge(edge.i, edge.j, weight);
    }

    NetworKit::APSP apsp(g);
    apsp.run();

    std::shared_ptr<std::ostream> out;

    if (param.outfile == "-") {
      out.reset(&std::cout, no_delete);
    } else {
      out =
        std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    for (uint64_t i = 1; i < h.attr.nodes; i++) {
      for (uint64_t j = 0; j < i; j++) {
        (*out) << apsp.getDistance(i, j) << ((j == i - 1) ? '\n' : '\t');
      }
    }

    rf.close();

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
