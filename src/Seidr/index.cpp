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
#include <index.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int
index(const std::vector<std::string>& args)
{
  logger log(std::cerr, "index");

  try {

    seidr_index_param_t param;

    po::options_description umbrella("Create index for SeidrFiles");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message");

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

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);

    param.infile = to_absolute(param.infile);
    param.outfile = param.infile + ".sfi";

    assert_exists(param.infile);
    assert_can_read(param.infile);
    assert_dir_is_writeable(dirname(param.outfile));

    log(LOG_INFO) << "Building index\n";

    SeidrFileIndex sfi;
    sfi.build(param.infile);

    log(LOG_INFO) << "Writing index to " << param.outfile << '\n';

    SeidrFile otf(param.outfile.c_str());
    otf.open("w");
    sfi.serialize(otf);
    otf.close();

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