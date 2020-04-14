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
#include <neighbours.h>
#include <parallel_control.h>
// External
#include <string>
#include <vector>

// Parallel includes
#if defined(SEIDR_PSTL)
#include <pstl/algorithm>
#include <pstl/execution>
#include <tbb/task_scheduler_init.h>
#else
#include <algorithm>
#endif

#include <armadillo>
#include <boost/numeric/conversion/cast.hpp>
#include <set>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int
neighbours(int argc, char* argv[])
{

  logger log(std::cerr, "neighbours");

  const char* args[argc - 1];
  std::string pr(argv[0]);
  pr += " neighbours";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++)
    args[i - 1] = argv[i];
  argc--;

  seidr_neighbours_param_t param;

  po::options_description umbrella(
    "Get the top N first-degree neighbours for each node or a"
    " list of nodes");

  po::options_description opt("Common Options");
  opt.add_options()("force,f",
                    po::bool_switch(&param.force)->default_value(false),
                    "Force overwrite output file if it exists")(
    "strict,s",
    po::bool_switch(&param.strict)->default_value(false),
    "Fail on all issues instead of warning")(
    "case-insensitive,I",
    po::bool_switch(&param.case_insensitive)->default_value(false),
    "Search case insensitive nodes")("help,h", "Show this help message")(
    "outfile,o",
    po::value<std::string>(&param.outfile)->default_value("-"),
    "Output file name ['-' for stdout]");

  po::options_description ompopt("OpenMP Options");
  ompopt.add_options()(
    "threads,O",
    po::value<int>(&param.nthreads)->default_value(GET_MAX_PSTL_THREADS()),
    "Number of OpenMP threads for parallel sorting");

  po::options_description nopt("Neighbours Options");
  nopt.add_options()(
    "index,i",
    po::value<uint32_t>(&param.tpos)->default_value(0, "last score"))(
    "neighbours,n",
    po::value<uint16_t>(&param.n)->default_value(10),
    "Number of top first-degree neighbours to return")(
    "genes,g", po::value<std::string>(&param.nodelist), "Gene names to search")(
    "rank,r",
    po::bool_switch(&param.trank)->default_value(false),
    "Use rank instead of score")(
    "unique,u",
    po::bool_switch(&param.unique)->default_value(false),
    "Print only unique edges");

  po::options_description req("Required [can be positional]");
  req.add_options()("in-file",
                    po::value<std::string>(&param.infile)->required(),
                    "Input SeidrFile");

  umbrella.add(req).add(nopt).add(ompopt).add(opt);

  po::positional_options_description p;
  p.add("in-file", 1);

  po::variables_map vm;
  po::store(
    po::command_line_parser(argc, args).options(umbrella).positional(p).run(),
    vm);

  if (vm.count("help") || argc == 1) {
    std::cerr << umbrella << '\n';
    return EINVAL;
  }

  try {
    po::notify(vm);
  } catch (std::exception& e) {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  try {
    set_pstl_threads(param.nthreads);

    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);

    param.indexfile = param.infile + ".sfi";
    assert_exists(param.indexfile);
    assert_can_read(param.indexfile);

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

  SeidrFile f(param.infile.c_str());
  f.open("r");
  SeidrFileHeader h;
  h.unserialize(f);

  make_tpos(param.tpos, h);

  std::shared_ptr<std::ostream> out;
  if (param.outfile == "-")
    out = std::shared_ptr<std::ostream>(&std::cout, [](void*) {});
  else
    out =
      std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));

  log(LOG_INFO) << "Reading indexfile...\n";
  SeidrFile sfi(param.indexfile.c_str());
  sfi.open("r");
  SeidrFileIndex index(param.case_insensitive, param.strict);
  index.unserialize(sfi, h);
  sfi.close();

  try {
    if (vm.count("genes")) {
      std::vector<std::string> genes = tokenize_delim(param.nodelist, ",");
      std::vector<SeidrFileEdge> edges_nu;
      std::set<SeidrFileEdge, sfe_score_sort_class> edges_u;
      for (auto& gene : genes) {
        std::vector<offset_t> offsets = index.get_offset_node(gene);
        if (offsets.size() == 0) {
          log(LOG_WARN) << "No edges found for " << gene << '\n';
        }
        std::vector<SeidrFileEdge> edges;
        for (auto& offset : offsets) {
          if (offset.o == -1) {
            log(LOG_WARN) << offset.err << '\n';
          } else {
            f.seek(offset.o);
            SeidrFileEdge e;
            e.unserialize(f, h);
            e.index.i = offset.i;
            e.index.j = offset.j;
            edges.push_back(e);
          }
        }
        if (param.trank) {
          SORTWCOMP(edges.begin(), edges.end(), sfe_score_sort);
        } else {
          SORTWCOMP(edges.begin(), edges.end(), sfe_rank_sort);
        }
        for (uint16_t i = 0; i < param.n; i++) {
          if (param.unique && i < edges.size()) {
            edges_u.insert(edges[i]);
          } else if (i < edges.size()) {
            edges_nu.push_back(edges[i]);
          }
        }
      }
      if (param.unique) {
        for (auto& e : edges_u) {
          e.print(h);
        }
      } else {
        for (auto& e : edges_nu) {
          e.print(h);
        }
      }
    } else {
      log(LOG_INFO) << "Creating score matrix...\n";
      arma::mat gm = read_network_arma(h,
                                       f,
                                       -std::numeric_limits<double>::infinity(),
                                       param.tpos,
                                       param.trank);
      std::vector<offset_t> offsets_nu;
      std::set<offset_t> offsets_u;
      for (arma::uword i = 0; i < gm.n_rows; i++) {
        arma::vec v = gm.col(i);

        arma::uvec sort_index;
        if (param.trank) {
          sort_index = arma::sort_index(v);
        } else {
          sort_index = arma::sort_index(v, "descend");
        }

        for (uint16_t j = 0; j < param.n; j++) {
          uint32_t xi = i;
          uint32_t xj = sort_index[j];
          offset_t off = index.get_offset_pair(h.nodes[xi], h.nodes[xj]);
          if (off.o == -1) {
            log(LOG_WARN) << off.err << '\n';
          } else if (param.unique) {
            offsets_u.insert(off);
          } else {
            offsets_nu.push_back(off);
          }
        }
      }
      if (param.unique) {
        for (auto& offset : offsets_u) {
          f.seek(offset.o);
          SeidrFileEdge e;
          e.unserialize(f, h);
          e.index.i = offset.i;
          e.index.j = offset.j;
          e.print(*out, h);
        }
      } else {
        for (auto& offset : offsets_nu) {
          f.seek(offset.o);
          SeidrFileEdge e;
          e.unserialize(f, h);
          e.index.i = offset.i;
          e.index.j = offset.j;
          e.print(*out, h);
        }
      }
    }
  } catch (std::exception& e) {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }
  return 0;
}