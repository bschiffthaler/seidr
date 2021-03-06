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
#include <stats.h>

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#undef DEBUG
#undef INFO

#include <boost/program_options.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/centrality/Closeness.hpp>
#include <networkit/centrality/DegreeCentrality.hpp>
#include <networkit/centrality/EigenvectorCentrality.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include <networkit/centrality/KatzCentrality.hpp>
#include <networkit/centrality/LaplacianCentrality.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/centrality/SpanningEdgeCentrality.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>

namespace po = boost::program_options;

bool
in_string(const std::string& query, const std::string& target)
{
  auto p = target.find(query);
  return p != std::string::npos;
}

int
stats(const std::vector<std::string>& args)
{
  logger log(std::cerr, "stats");

  // Turn off NetworKit's logging
  Aux::Log::setLogLevel("QUIET");

  try {

    seidr_stats_param_t param;

    po::options_description umbrella("Calculate network centrality statistics");

    po::options_description opt("Common Options");
    opt.add_options()("help,h", "Show this help message")(
      "tempdir,T",
      po::value<std::string>(&param.tempdir)->default_value("", "auto"),
      "Directory to store temporary data");

    po::options_description sopt("Stats Options");
    sopt.add_options()(
      "index,i",
      po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
      "Index of score to use")(
      "nsamples,n",
      po::value<uint32_t>(&param.nsamples)->default_value(0),
      "Use N samples for approximations")(
      "metrics,m",
      po::value<std::string>(&param.metrics)
        ->default_value("PR,CLO,BTW,STR,EV,KTZ,LPC,SEC,EBC"),
      "String describing metrics to calculate")(
      "directed,d",
      po::bool_switch(&param.directed)->default_value(false),
      "(Experimental) Use directionality information.")(
      "eigenvector-tol",
      po::value<double>(&param.ev_tol)->default_value(SEIDR_STATS_DEF_EV_TOL),
      "Eigenvector centrality convergence tolerance")(
      "exact,e",
      po::bool_switch(&param.exact)->default_value(false),
      "Calculate exact stats.")(
      "weight-is-distance",
      po::bool_switch(&param.w_is_dist)->default_value(false),
      "Edge weight represents a distance (rather than similarity).")(
      "weight-rank",
      po::bool_switch(&param.trank)->default_value(false),
      "Set weight value to rank rather than score");

    po::options_description req("Required");
    req.add_options()("in-file",
                      po::value<std::string>(&param.infile)->required(),
                      "Input SeidrFile [can be positional]");

    umbrella.add(req).add(sopt).add(opt);

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
    for (const auto& metric : tokenize_delim(param.metrics, ",")) {
      assert_arg_constraint<std::string>(
        { "PR", "CLO", "BTW", "STR", "EV", "KTZ", "SEC", "EBC", "LPC" },
        metric);
    }

    SeidrFile rf(param.infile.c_str());
    rf.open("r");

    SeidrFileHeader h;
    h.unserialize(rf);

    // if (h.attr.pagerank_calc == 1 && in_string("PR", param.metrics)) {
    //   log(LOG_WARN) << "Found previous PageRank calculation. Deleting...\n";
    //   h.attr.pagerank_calc = 0;
    //   h.pagerank.clear();
    // }
    // if (h.attr.closeness_calc == 1 && in_string("CLO", param.metrics)) {
    //   log(LOG_WARN) << "Found previous closeness calculation. Deleting...\n";
    //   h.attr.closeness_calc = 0;
    //   h.closeness.clear();
    // }
    // if (h.attr.betweenness_calc == 1 && in_string("BTW", param.metrics)) {
    //   log(LOG_WARN) << "Found previous betweenness calculation.
    //   Deleting...\n"; h.attr.betweenness_calc = 0; h.betweenness.clear();
    // }
    // if (h.attr.strength_calc == 1 && in_string("STR", param.metrics)) {
    //   log(LOG_WARN) << "Found previous strength calculation. Deleting...\n";
    //   h.attr.strength_calc = 0;
    //   h.strength.clear();
    // }
    // if (h.attr.eigenvector_calc == 1 && in_string("EV", param.metrics)) {
    //   log(LOG_WARN) << "Found previous eigenvector calculation.
    //   Deleting...\n"; h.attr.eigenvector_calc = 0; h.eigenvector.clear();
    // }
    // if (h.attr.katz_calc == 1 && in_string("KTZ", param.metrics)) {
    //   log(LOG_WARN) << "Found previous katz calculation. Deleting...\n";
    //   h.attr.katz_calc = 0;
    //   h.katz.clear();
    // }

    std::vector<std::string> centrality_names;
    std::vector<double> centrality_data;

    if (!param.exact && param.nsamples == 0) {
      param.nsamples = boost::numeric_cast<uint64_t>(
        boost::numeric_cast<double>(h.attr.nodes) * SEIDR_STATS_DEF_SAMPLES);
      log(LOG_INFO) << "Using " << param.nsamples
                    << " as the number of samples\n";
    }
    // Workaround for GitHub issue #9 (crash because approximate data structure
    // is initialized even though we are not using it)
    if (param.nsamples == 0) {
      param.nsamples = h.attr.nodes;
    }

    make_tpos(param.tpos, h);

    log(LOG_INFO) << "Reading graph\n";
    std::vector<SeidrFileEdge> ev = read_network(
      h, rf, -std::numeric_limits<double>::infinity(), param.tpos, false);

    log(LOG_INFO) << "Creating Networkit graph\n";
    NetworKit::Graph g(h.attr.nodes, true, param.directed);
    uint64_t ud = 0;
    uint64_t ab = 0;
    uint64_t ba = 0;
    for (auto& e : ev) {
      double weight =
        param.trank ? e.scores[param.tpos].r : e.scores[param.tpos].s;
      // We do higher=better calculations first, so we need to inverse weight
      // if it represents a distance
      weight = param.w_is_dist ? 1 / weight : weight;
      if (param.directed) {
        if (EDGE_IS_DIRECT(e.attr.flag)) {
          if (EDGE_IS_AB(e.attr.flag)) {
            g.addEdge(e.index.i, e.index.j, weight);
            ab++;
          } else {
            g.addEdge(e.index.j, e.index.i, weight);
            ba++;
          }
        } else {
          // No direction so we add both i->j and j->i
          g.addEdge(e.index.i, e.index.j, weight);
          g.addEdge(e.index.j, e.index.i, weight);
          ab++;
          ba++;
        }
      } else {
        g.addEdge(e.index.i, e.index.j, weight);
        ud++;
      }
    }

    g.indexEdges();

    log(LOG_INFO) << "Have graph with " << g.numberOfNodes() << " nodes and "
                  << g.numberOfEdges() << " edges\n";
    if (param.directed) {
      log(LOG_INFO) << ab << " edges in direction A->B, " << ba << " edges "
                    << "in direction A<-B\n";
    }

    rf.close();

    log(LOG_INFO) << "Calculating connected components\n";
    uint64_t comp = 0;
    if (param.directed) {
      NetworKit::StronglyConnectedComponents ccn(g);
      ccn.run();
      comp = ccn.numberOfComponents();
    } else {
      NetworKit::ConnectedComponents ccn(g);
      ccn.run();
      comp = ccn.numberOfComponents();
    }

    log(LOG_INFO) << "Have " << comp << " components\n";

    if (in_string("PR", param.metrics)) {
      log(LOG_INFO) << "Calculating PageRank\n";
      NetworKit::PageRank prv(g);
      prv.run();
      // h.attr.pagerank_calc = 1;
      centrality_names.emplace_back("PageRank");
      for (uint32_t i = 0; i < h.attr.nodes; i++) {
        centrality_data.push_back(prv.score(i));
      }
    }

    if (in_string("STR", param.metrics)) {
      log(LOG_INFO) << "Calculating Node Strength\n";
      // h.attr.strength_calc = 1;
      if (param.directed) {
        centrality_names.emplace_back("OutStrength");
      } else {
        centrality_names.emplace_back("Strength");
      }
      
      std::vector<double> indeg;
      for (uint64_t i = 0; i < h.attr.nodes; i++) {
        double xdeg = g.weightedDegree(i);
        double xindeg = g.weightedDegreeIn(i);
        centrality_data.push_back(xdeg);
        indeg.push_back(xindeg);
      }
      if (param.directed) {
        centrality_names.emplace_back("InStrength");
        for (uint64_t i = 0; i < h.attr.nodes; i++) {
          centrality_data.push_back(indeg[i]);
        }
      }
    }

    if (in_string("EV", param.metrics)) {
      log(LOG_INFO) << "Calculating Eigenvector centrality\n";
      NetworKit::EigenvectorCentrality evv(g, param.ev_tol);
      evv.run();
      // h.attr.eigenvector_calc = 1;
      centrality_names.emplace_back("Eigenvector");
      for (uint32_t i = 0; i < h.attr.nodes; i++) {
        centrality_data.push_back(evv.score(i));
      }
    }

    if (in_string("KTZ", param.metrics)) {
      log(LOG_INFO) << "Calculating Katz centrality\n";
      NetworKit::KatzCentrality kav(g);
      kav.run();
      // h.attr.katz_calc = 1;
      centrality_names.emplace_back("Katz");
      for (uint32_t i = 0; i < h.attr.nodes; i++) {
        centrality_data.push_back(kav.score(i));
      }
    }

    if (in_string("LPC", param.metrics)) {
      log(LOG_INFO) << "Calculating Laplacian centrality\n";
      NetworKit::LaplacianCentrality lav(g);
      lav.run();
      // h.attr.katz_calc = 1;
      centrality_names.emplace_back("Laplacian");
      for (uint32_t i = 0; i < h.attr.nodes; i++) {
        centrality_data.push_back(lav.score(i));
      }
    }

    // If weight is a distance it was inversed before and we need to now
    // undo that. If it's a similarity we need to turn it into a distance.
    // In any case, we inverse again
    for (const auto e : g.edgeWeightRange()) {
      g.setWeight(e.u, e.v, 1.0 / e.weight);
    }

    if (in_string("CLO", param.metrics)) {
      log(LOG_INFO) << "Calculating " << (param.exact ? "exact" : "approximate")
                    << " Closeness\n";

      if (comp > 1) {
        log(LOG_WARN)
          << "Graph is disconnected, switching to exact variant...\n";
      }

      NetworKit::ApproxCloseness clv(g, param.nsamples);
      NetworKit::Closeness clve(
        g, false, NetworKit::ClosenessVariant::generalized);
      if (param.exact || comp > 1) {
        clve.run();
      } else {
        clv.run();
      }
      for (uint32_t i = 0; i < h.attr.nodes; i++) {
        centrality_data.push_back((param.exact || comp > 1) ? clve.score(i) : clv.score(i));
      }
      centrality_names.emplace_back("Closeness");
    }

    std::vector<double> ebv;
    if (in_string("BTW", param.metrics)) {
      if (!param.exact && in_string("EBC", param.metrics)) {
        log(LOG_WARN) << "Edge betweenness requires exact statistics. Forcing "
                         "exact betweenness\n";
      } else {
        log(LOG_INFO) << "Calculating "
                      << (param.exact ? "exact" : "approximate")
                      << " Betweenness\n";
      }
      NetworKit::EstimateBetweenness btv(g, param.nsamples, false, true);
      NetworKit::Betweenness btve(g, false, in_string("EBC", param.metrics));
      if (param.exact || in_string("EBC", param.metrics)) {
        btve.run();
        if (in_string("EBC", param.metrics)) {
          ebv = btve.edgeScores();
        }
      } else {
        btv.run();
      }
      for (uint32_t i = 0; i < h.attr.nodes; i++) {
        centrality_data.push_back(
          (param.exact || in_string("EBC", param.metrics)) ? btve.score(i)
                                                           : btv.score(i));
      }
      centrality_names.emplace_back("Betweenness");
      // h.attr.betweenness_calc = 1;
    }

    std::vector<double> spv;
    if (in_string("SEC", param.metrics)) {
      log(LOG_INFO) << "Calculating Spanning Edge centrality\n";

      NetworKit::SpanningEdgeCentrality spc(g);
      if (param.exact) {
        spc.run();
      } else {
        spc.runParallelApproximation();
      }
      spv = spc.scores();
    }

    log(LOG_INFO) << "Writing results\n";

    bool sec_found = false;
    bool ebc_found = false;
    size_t sec_i = 0;
    size_t ebc_i = 0;

    if (in_string("SEC", param.metrics)) {
      auto sec_ptr = std::find(h.supp.begin(), h.supp.end(), "SEC");
      if (sec_ptr == h.supp.end()) {
        h.attr.nsupp += 1;
        h.attr.nsupp_flt += 1;
        h.supp.emplace_back("SEC");
      } else {
        sec_found = true;
        sec_i = std::distance(h.supp.begin(), sec_ptr);
        sec_i -= h.attr.nsupp_str;
        sec_i -= h.attr.nsupp_int;
      }
    }

    if (in_string("EBC", param.metrics)) {
      auto ebc_ptr = std::find(h.supp.begin(), h.supp.end(), "EBC");
      if (ebc_ptr == h.supp.end()) {
        h.attr.nsupp += 1;
        h.attr.nsupp_flt += 1;
        h.supp.emplace_back("EBC");
      } else {
        ebc_found = true;
        ebc_i = std::distance(h.supp.begin(), ebc_ptr);
        ebc_i -= h.attr.nsupp_str;
        ebc_i -= h.attr.nsupp_int;
      }
    }

    param.tempfile = tempfile(param.tempdir);
    SeidrFile tf(param.tempfile.c_str());
    tf.open("w");

    h.centrality_data = centrality_data;
    h.centrality_names = centrality_names;
    h.attr.pagerank_calc = centrality_names.size();
    // We update header version or else bad stuff happens
    h.version_from_char(_XSTR(VERSION));
    h.serialize(tf);

    for (uint64_t i = 0; i < ev.size(); i++) {
      SeidrFileEdge x = ev[i];
      if (in_string("SEC", param.metrics)) {
        if (sec_found) {
          x.supp_flt[sec_i] = spv[i];
        } else {
          x.supp_flt.push_back(spv[i]);
        }
      }
      if (in_string("EBC", param.metrics)) {
        if (ebc_found) {
          x.supp_flt[ebc_i] = ebv[i];
        } else {
          x.supp_flt.push_back(ebv[i]);
        }
      }
      x.serialize(tf, h);
    }

    tf.close();

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
