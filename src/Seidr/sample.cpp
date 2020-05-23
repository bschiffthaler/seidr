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
#include <cmath>
#include <common.h>
#include <fs.h>
#include <sample.h>
// External
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <bs/simple_sample.h>
#include <bs/vitter_d.h>
#include <random>

namespace po = boost::program_options;

uint64_t
prop_to_nedges(double prop, uint64_t n)
{
  auto nd = boost::numeric_cast<double>(n) * prop;
  auto ret = boost::numeric_cast<uint64_t>(round(nd));
  return ret;
}

void
sample_wo_repl(const seidr_sample_param_t& param,
               const po::variables_map& vm,
               std::ostream& out)
{
  SeidrFile sf(param.el_file.c_str());
  sf.open("r");
  SeidrFileHeader header;
  header.unserialize(sf);
  sf.close();
  sf.open("r");
  uint64_t nedges = 0;
  if (vm.count("nedges") != 0) {
    if (param.nedges > header.attr.edges) {
      throw std::runtime_error("Can't sample " + std::to_string(param.nedges) +
                               " edges when the network only contains " +
                               std::to_string(header.attr.edges) + " edges");
    }
    nedges = param.nedges;
  } else {
    nedges = prop_to_nedges(param.prop, header.attr.edges);
  }

  uint64_t ctr = 0;

  BS::vitter_d vd(header.attr.edges, nedges);

  uint64_t cur = vd.next();

  sf.each_edge_exit_early([&](SeidrFileEdge& e, SeidrFileHeader& h) {
    if (cur == ctr) {
      e.print(out, h);
      if (vd.end()) {
        return true;
      }
      cur = vd.next();
    }
    ctr++;
    return false;
  });

  sf.close();
}

void
sample_w_repl(const seidr_sample_param_t& param,
              const po::variables_map& vm,
              std::ostream& out)
{
  SeidrFile sf(param.el_file.c_str());
  sf.open("r");
  SeidrFileHeader header;
  header.unserialize(sf);
  sf.close();
  sf.open("r");
  uint64_t nedges = 0;
  if (vm.count("nedges") != 0) {
    nedges = param.nedges;
  } else {
    nedges = prop_to_nedges(param.prop, header.attr.edges);
  }

  uint64_t ctr = 0;

  BS::simple_sample<uint64_t> ss(header.attr.edges, nedges);

  uint64_t cur = ss.next();

  sf.each_edge_exit_early([&](SeidrFileEdge& e, SeidrFileHeader& h) {
    while (cur == ctr) {
      e.print(out, h);
      if (ss.end()) {
        return true;
      }
      cur = ss.next();
    }
    ctr++;
    return false;
  });

  sf.close();
}

int
sample(const std::vector<std::string>& args)
{

  logger log(std::cerr, "sample");

  try {

    // Variables used by the function
    seidr_sample_param_t param;

    po::options_description umbrella;

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message")(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("-"),
      "Output file name ['-' for stdout]");

    po::options_description sopt("Sample options");
    sopt.add_options()(
      "replacement,r",
      po::bool_switch(&param.replacement)->default_value(false),
      "Sample with replacement")(
      "precision,p",
      po::value<uint16_t>(&param.prec)->default_value(SEIDR_SAMPLE_DEF_PREC),
      "Number of significant digits to report")(
      "nedges,n",
      po::value<uint64_t>(&param.nedges),
      "Number of edges to sample")("fraction,F",
                                   po::value<double>(&param.prop),
                                   "Fraction of edges to sample");

    po::options_description req("Required");
    req.add_options()("in-file",
                      po::value<std::string>(&param.el_file)->required(),
                      "Input SeidrFile");

    umbrella.add(req).add(sopt).add(opt);

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

    // Check a bunch of sanity things
    assert_exists(param.el_file);
    assert_can_read(param.el_file);

    if (vm.count("nedges") != 0 && vm.count("fraction") != 0) {
      throw std::runtime_error(
        "Supply either --fraction or --nedges, not both");
    }

    if ((vm.count("nedges") == 0) && (vm.count("fraction") == 0)) {
      throw std::runtime_error("Supply either --fraction or --nedges");
    }

    if (vm.count("fraction") != 0) {
      if (!param.replacement && (param.prop < 0 || param.prop > 1)) {
        throw std::runtime_error("Fraction must be in [0, 1] if sampling "
                                 "without replacement");
      }
      if (param.replacement && param.prop < 0) {
        throw std::runtime_error("Fraction must be > 0");
      }
    }

    std::shared_ptr<std::ostream> out;
    if (param.outfile == "-") {
      out = std::shared_ptr<std::ostream>(&std::cout, no_delete);
    }
    else {
      out =
        std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    // Start main function
    if (!param.replacement) {
      sample_wo_repl(param, vm, *out);
    } else {
      sample_w_repl(param, vm, *out);
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
