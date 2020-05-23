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
#include <common.h>
#include <fs.h>
#include <reheader.h>

#include <boost/program_options.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace po = boost::program_options;

int
reheader(const std::vector<std::string>& args)
{
  logger log(std::cerr, "reheader");

  try {

    seidr_reheader_param_t param;

    po::options_description umbrella(
      "Modify SeidrFile headers.\n"
      "Currently only drops disconnected nodes and resets\n"
      "stats.");

    po::options_description opt("Common Options");
    opt.add_options()("help,h", "Show this help message")(
      "tempdir,T",
      po::value<std::string>(&param.tempdir)->default_value("", "auto"),
      "Directory to store temporary data");

    po::options_description req("Required [can be positional]");
    req.add_options()("in-file",
                      po::value<std::string>(&param.infile)->required(),
                      "Input SeidrFile");

    umbrella.add(req).add(opt);

    po::positional_options_description p;
    p.add("in-file", 1);

    po::variables_map vm;
    po::store(
      po::command_line_parser(args).options(umbrella).positional(p).run(),
      vm);

    if (vm.count("help") != 0) {
      std::cerr << umbrella << '\n';
      return EINVAL;
    }

    po::notify(vm);

    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);
    assert_dir_is_writeable(param.tempdir);

    param.tempfile = tempfile(param.tempdir);

    SeidrFile rf(param.infile.c_str());
    rf.open("r");

    // Fill map with old index data
    std::map<uint32_t, std::string> old_index;
    rf.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h) {
      auto hit = old_index.find(e.index.i);
      if (hit == old_index.end()) {
        old_index[e.index.i] = h.nodes[e.index.i];
      }
      hit = old_index.find(e.index.j);
      if (hit == old_index.end()) {
        old_index[e.index.j] = h.nodes[e.index.j];
      }
    });

    rf.close();

    // Create new vector of nodes for the header
    std::vector<std::string> new_nodes;
    new_nodes.reserve(old_index.size());
    for (const auto& e : old_index) {
      new_nodes.push_back(e.second);
    }

    // Fill a map with new index data
    std::map<std::string, uint32_t> new_index;
    for (uint32_t i = 0; i < new_nodes.size(); i++) {
      new_index[new_nodes[i]] = i;
    }

    SeidrFile of(param.tempfile.c_str());
    of.open("w");

    SeidrFile inf(param.infile.c_str());
    inf.open("r");

    SeidrFileHeader h;
    h.unserialize(inf);

    auto old_nodes = h.nodes;

    // Fix nodes
    h.attr.nodes = new_nodes.size();
    h.nodes = new_nodes;
    h.serialize(of);

    inf.seek(0);
    // NOLINTNEXTLINE: SeidrFileHeader& must be declared even if not used as a parameter
    inf.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h2) { 
      e.index.i = new_index[old_index.at(e.index.i)];
      e.index.j = new_index[old_index.at(e.index.j)];
      e.serialize(of, h); // Use NEW header
    });

    of.close();
    inf.close();

    rename(param.tempfile, param.infile);

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