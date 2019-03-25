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
#include <common.h>
#include <resolve.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <cerrno>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

void resolve_im(SeidrFileHeader& h, std::string& fname, std::ostream& out)
{

  std::ifstream ifs(fname.c_str(), std::ios::in);
  std::string line;
  uint32_t sc_max = 0;
  while (std::getline(ifs, line))
  {
    if (line.at(0) == '#') continue;
    std::string field;
    std::stringstream ss(line);
    ss >> field;
    uint32_t sc = 0;
    for (char& c : field)
      if (c == ':')
        sc++;
    sc_max = sc > sc_max ? sc : sc_max;
  }

  ifs.clear();
  ifs.seekg(0, std::ios_base::beg);

  while (std::getline(ifs, line))
  {
    if (line.at(0) == '#') continue;
    std::string field;
    std::stringstream ss(line);
    uint32_t ctr = 0;
    while (std::getline(ss, field, ' '))
    {
      if (ctr == 0)
      {
        out << field << '\t';
        uint32_t sc = 0;
        for (char c : field)
          if (c == ':')
          {
            out << '\t';
            sc++;
          }
          else
            out << c;
        while (sc < sc_max)
        {
          out << '\t' << "NA";
          sc++;
        }
        out << '\t';
        ctr++;
      }
      else if (ctr == 3)
      {
        auto index = std::stoul(field);
        out << h.nodes[index] << '\n';
        ctr = 0;
      }
      else
      {
        out << field << '\t';
        ctr++;
      }
    }
  }
}

int resolve(int argc, char * argv[])
{

  logger log(std::cerr, "resolve");

  // We ignore the first argument, the function name
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " resolve";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  seidr_resolve_param_t param;

  po::options_description
  umbrella("Resolve node indices in text file to node names.");

  po::options_description opt("Common Options");
  opt.add_options()
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("help,h", "Show this help message")
  ("outfile,o",
   po::value<std::string>(&param.outfile)->default_value("-"),
   "Output file name ['-' for stdout]");

  po::options_description ropt("Resolve Options");
  ropt.add_options()
  ("seidr-file,s",
   po::value<std::string>(&param.sf),
   "Seidr file which should be used to resolve input")
  ("format,F",
   po::value<std::string>(&param.format)->default_value("infomap"),
   "File format to resolve");

  po::options_description req("Required [can be positional]");
  req.add_options()
  ("in-file", po::value<std::string>(&param.infile)->required(),
   "Input SeidrFile");

  umbrella.add(req).add(ropt).add(opt);

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

  try
  {
    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);

    if (param.outfile != "-")
    {
      param.outfile = to_absolute(param.outfile);
      if (! param.force)
      {
        assert_no_overwrite(param.outfile);
      }
      assert_dir_is_writeable(dirname(param.outfile));
    }
    assert_arg_constraint<std::string>({"infomap"}, param.format);
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  std::shared_ptr<std::ostream> out;
  if (param.outfile == "-")
    out = std::shared_ptr < std::ostream > (&std::cout, [](void*) {});
  else
    out = std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));

  if (param.format == "infomap")
  {
    SeidrFile fin(param.sf.c_str());
    fin.open("r");
    SeidrFileHeader h;
    h.unserialize(fin);
    resolve_im(h, param.infile, *out);
  }

  return 0;
}
