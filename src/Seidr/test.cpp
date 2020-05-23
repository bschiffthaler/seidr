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
#include <test.h>
// External
#include <cerrno>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

int
test(const std::vector<std::string>& args)
{

  // logger log(std::cerr, "test");

  // // Variables used by the function
  // std::string el_file;

  // try
  // {
  //   el_file = to_canonical(el_file);
  // }
  // catch (std::runtime_error& e)
  // {
  //   log(LOG_ERR) << e.what() << '\n';
  //   return errno;
  // }

  // SeidrFile rf(el_file.c_str());
  // rf.open("r");

  // rf.each_edge([](SeidrFileEdge& e, SeidrFileHeader& h){
  //   e.print(h);
  // });

  // rf.close();

  return 0;
}
