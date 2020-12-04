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
#include <top.h>
// External
#include <boost/program_options.hpp>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace po = boost::program_options;

typedef std::pair<SeidrFileEdge, double> expl_score_t;

// Priority queue should be sorted inverse (greater scores last) to make
// popping off the minimum element easy
class comp_greater
{
public:
  bool operator()(const expl_score_t& lhs, const expl_score_t& rhs)
  {
    return lhs.second > rhs.second;
  }
};

// Need to derive from STL priority_queue to access the underlying vector
// for convenient iteration
class priority_queue
  : public std::
      priority_queue<expl_score_t, std::vector<expl_score_t>, comp_greater>
{
public:
  std::vector<expl_score_t>& vec() { return c; }
};

int
top(const std::vector<std::string>& args)
{
  logger log(std::cerr, "top");

  try {

    seidr_top_param_t param;

    po::options_description umbrella("Get edges with highest scores in a "
                                     "memory friendly way.");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("-"),
      "Output file name ['-' for stdout]");

    po::options_description topt("Top Options");
    topt.add_options()(
      "number,n",
      po::value<uint64_t>(&param.ntop)->default_value(SEIDR_TOP_DEF_NTOP),
      "The number of highest scoring edges to return")(
      "index,i",
      po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
      "Score column to use as edge weights");

    po::options_description req("Required [can be positional]");
    req.add_options()("in-file",
                      po::value<std::string>(&param.infile)->required(),
                      "Input SeidrFile");

    umbrella.add(req).add(topt).add(opt);

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

    SeidrFile rf(param.infile.c_str());
    rf.open("r");
    SeidrFileHeader h;
    h.unserialize(rf);

    rf.seek(0);

    make_tpos(param.tpos, h);

    double min = -std::numeric_limits<double>::infinity();
    double max = min;

    std::shared_ptr<std::ostream> out;
    if (param.outfile == "-") {
      out = std::shared_ptr<std::ostream>(&std::cout, no_delete);
    } else {
      out =
        std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    priority_queue pq;

    // NOLINTNEXTLINE: SeidrFileHeader& needs to be declared even if unused
    rf.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h) {
      double v = e.scores[param.tpos].s;
      if (pq.size() < param.ntop) {
        pq.push(expl_score_t(e, v));
        if (v > max) {
          max = v;
        }
        min = pq.top().second;
      } else if (v > min || almost_equal(v, min)) {
        pq.push(expl_score_t(e, v));
        if (v > max) {
          max = v;
        }
        min = pq.top().second;
      }
      if (pq.size() > param.ntop) {
        pq.pop();
        min = pq.top().second;
      }
    });

    while (!pq.empty()) {
      pq.top().first.print(*out, h);
      pq.pop();
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
