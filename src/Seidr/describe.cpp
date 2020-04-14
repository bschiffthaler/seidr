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
#include <describe.h>
#include <fs.h>
// External
#include <armadillo>
#include <boost/program_options.hpp>
#include <bs/describe.h>
#include <bs/histogram.h>
#include <cerrno>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace po = boost::program_options;

void
get_desc_stats(seidr_param_describe_t& params)
{

  std::shared_ptr<std::ostream> out;
  if (params.outfile == "-") {
    out.reset(&std::cout, [](...) {});
  } else {
    out.reset(new std::ofstream(params.outfile.c_str()));
  }

  SeidrFile sf(params.el_file.c_str());
  sf.open("r");

  SeidrFileHeader header;
  header.unserialize(sf);
  sf.seek(0);

  std::vector<BS::desc_stats> stats;
  for (uint64_t i = 0; i < header.attr.nalgs; i++) {
    stats.push_back(BS::desc_stats());
  }

  sf.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h) {
    for (uint16_t i = 0; i < h.attr.nalgs; i++) {
      if (std::isfinite(e.scores[i].s)) {
        stats[i].add(e.scores[i].s);
      }
    }
  });

  sf.close();

  for (uint64_t i = 0; i < header.attr.nalgs; i++) {
    (*out) << i << '\t' << header.algs[i] << "\tmin\t" << stats[i].min()
           << '\n';
    (*out) << i << '\t' << header.algs[i] << "\tmax\t" << stats[i].max()
           << '\n';
    (*out) << i << '\t' << header.algs[i] << "\tmean\t" << stats[i].mean()
           << '\n';
    (*out) << i << '\t' << header.algs[i] << "\tn\t" << stats[i].size() << '\n';

    std::string quant;
    std::stringstream qs(params.quantiles);
    while (std::getline(qs, quant, ',')) {
      (*out) << i << '\t' << header.algs[i] << "\tQ" << quant << '\t'
             << stats[i].quantile(std::stod(quant)) << '\n';
    }

    BS::histogram hist(stats[i], params.bins);
    std::stringstream ss;
    if (params.parseable) {
      hist.print_tsv(ss);
    } else {
      hist.print_vertical(ss);
    }

    std::string line;
    for (uint64_t j = 0; j < params.bins; j++) {
      std::getline(ss, line);
      (*out) << i << '\t' << header.algs[i] << "\tHIST\t" << line << '\n';
    }
  }
}

int
describe(int argc, char* argv[])
{

  LOG_INIT_CERR();

  // Variables used by the function
  seidr_param_describe_t param;

  // We ignore the first argument
  const char* args[argc - 1];
  std::string pr(argv[0]);
  pr += " describe";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++)
    args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella("Calculate summary statistics for each "
                                   "algorithm in a SeidrFile");

  po::options_description opt("Common Options");
  opt.add_options()("force,f",
                    po::bool_switch(&param.force)->default_value(false),
                    "Force overwrite output file if it exists")(
    "help,h", "Show this help message")(
    "outfile,o",
    po::value<std::string>(&param.outfile)->default_value("-"),
    "Output file name ['-' for stdout]");

  po::options_description dopt("Describe Options");
  dopt.add_options()("bins,b",
                     po::value<uint32_t>(&param.bins)->default_value(25),
                     "Number of histgram bins")(
    "parseable,p",
    po::bool_switch(&param.parseable)->default_value(false),
    "Print parseable histogram")(
    "quantiles,q",
    po::value<std::string>(&param.quantiles)
      ->default_value("0.05,0.25,0.5,0.75,0.95"),
    "A comma seaprated string of quantiles [0.05,0.25,0.5,0.75,0.95]");

  po::options_description req("Required Options [can be positional]");
  req.add_options()("in-file",
                    po::value<std::string>(&param.el_file)->required(),
                    "Input SeidrFile");

  umbrella.add(req).add(dopt).add(opt);

  po::positional_options_description p;
  p.add("in-file", 1);

  po::variables_map vm;
  po::store(
    po::command_line_parser(argc, args).options(umbrella).positional(p).run(),
    vm);

  if (vm.count("help") || argc == 1) {
    std::cerr << umbrella << '\n';
    return 22;
  }

  try {
    po::notify(vm);
    param.el_file = to_absolute(param.el_file);
    assert_exists(param.el_file);
    assert_can_read(param.el_file);
    if (param.outfile != "-") {
      param.outfile = to_absolute(param.outfile);
      if (!param.force) {
        assert_no_overwrite(param.outfile);
      }
      assert_dir_is_writeable(dirname(param.outfile));
    }
  } catch (std::runtime_error& except) {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  get_desc_stats(param);

  return 0;
}