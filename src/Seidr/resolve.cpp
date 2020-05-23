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
#include <resolve.h>
// External
#include <boost/program_options.hpp>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace po = boost::program_options;

void
resolve_im(SeidrFileHeader& h, const std::string& fname, std::ostream& out)
{

  std::ifstream ifs(fname.c_str(), std::ios::in);
  std::string line;
  uint32_t sc_max = 0;
  while (std::getline(ifs, line)) {
    if (line.empty() || line.at(0) == '#') {
      continue;
    }
    std::string field;
    std::stringstream ss(line);
    ss >> field;
    uint32_t sc = 0;
    for (char& c : field) {
      if (c == ':') {
        sc++;
      }
    }
    sc_max = sc > sc_max ? sc : sc_max;
  }

  ifs.clear();
  ifs.seekg(0, std::ios_base::beg);

  while (std::getline(ifs, line)) {
    if (line.empty() || line.at(0) == '#') {
      continue;
    }
    std::string field;
    std::stringstream ss(line);
    uint32_t ctr = 0;
    while (std::getline(ss, field, ' ')) {
      if (ctr == 0) {
        out << field << '\t';
        uint32_t sc = 0;
        for (char c : field) {
          if (c == ':') {
            out << '\t';
            sc++;
          } else {
            out << c;
          }
        }
        while (sc < sc_max) {
          out << '\t' << "NA";
          sc++;
        }
        out << '\t';
        ctr++;
      } else if (ctr == 3) {
        auto index = std::stoul(field);
        out << h.nodes[index] << '\n';
        ctr = 0;
      } else {
        out << field << '\t';
        ctr++;
      }
    }
  }
}

int
resolve(const std::vector<std::string>& args)
{

  logger log(std::cerr, "resolve");

  try {

    seidr_resolve_param_t param;

    po::options_description umbrella(
      "Resolve node indices in text file to node names.");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("-"),
      "Output file name ['-' for stdout]");

    po::options_description ropt("Resolve Options");
    ropt.add_options()("seidr-file,s",
                       po::value<std::string>(&param.sf),
                       "Seidr file which should be used to resolve input")(
      "format,F",
      po::value<std::string>(&param.format)->default_value("infomap"),
      "File format to resolve");

    po::options_description req("Required [can be positional]");
    req.add_options()("in-file",
                      po::value<std::string>(&param.infile)->required(),
                      "Input SeidrFile");

    umbrella.add(req).add(ropt).add(opt);

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
    assert_arg_constraint<std::string>({ "infomap" }, param.format);

    std::shared_ptr<std::ostream> out;
    if (param.outfile == "-") {
      out = std::shared_ptr<std::ostream>(&std::cout, no_delete);
    } else {
      out =
        std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    if (param.format == "infomap") {
      SeidrFile fin(param.sf.c_str());
      fin.open("r");
      SeidrFileHeader h;
      h.unserialize(fin);
      resolve_im(h, param.infile, *out);
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
