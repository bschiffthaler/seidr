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

#include <common.h>
#include <index.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int index(int argc, char ** argv)
{
  logger log(std::cerr, "index");

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " index";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  seidr_index_param_t param;

  po::options_description umbrella("Create index for SeidrFiles");

  po::options_description opt("Common Options");
  opt.add_options()
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("help,h", "Show this help message");

  po::options_description req("Required [can be positional]");
  req.add_options()
  ("in-file", po::value<std::string>(&param.infile)->required(),
   "Input SeidrFile");

  umbrella.add(req).add(opt);

  po::positional_options_description p;
  p.add("in-file", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, args).
            options(umbrella).positional(p).run(), vm);

  if (vm.count("help") || argc == 1)
  {
    std::cerr << umbrella << '\n';
    return 22;
  }

  try
  {
    po::notify(vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  param.infile = to_absolute(param.infile);
  param.outfile = param.infile + ".sfi";

  try
  {
    assert_exists(param.infile);
    assert_can_read(param.infile);
    assert_dir_is_writeable(dirname(param.outfile));
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

  log(LOG_INFO) << "Building index\n";

  SeidrFileIndex sfi;
  sfi.build(param.infile);

  log(LOG_INFO) << "Writing index to " << param.outfile << '\n';

  SeidrFile otf(param.outfile.c_str());
  otf.open("w");
  sfi.serialize(otf);
  otf.close();
  return 0;
}