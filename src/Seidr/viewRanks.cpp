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
#include <viewRanks.h>
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <cerrno>
#include <iostream>
#include <memory>
#include <vector>
#include <stack>
#include <string>
#include <set>
#include <sstream>
#include <fstream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

static const std::map<char, uint8_t> op_prec_char =
{
  {'<', 10}, {'>', 10}, {'=', 9},
  {'&', 5}, {'|', 4}, {'(', 0},
  {')', 0}, {'+', 12}, {'-', 12},
  {'/', 13}, {'*', 13}
};

// The filter is implemented as a stack machine, so the easiest way is to
// first conver the filter string to postfix notation
std::vector<std::string> infix_to_postfix(std::string str)
{
  std::stack<char> opstack;
  std::vector<std::string> postfix;
  boost::char_separator<char> sep(" ", "<>=&|()+-/*");
  tokenizer tokens(str, sep);
  for (auto tok = tokens.begin(); tok != tokens.end(); tok++)
  {
    char tk = (*tok)[0];
    //std::cout << tk << '\n';
    auto t = op_prec_char.find(tk);
    if (t == op_prec_char.end())
    {
      //std::cout << "Push PF\n";
      postfix.push_back(*tok);
    }
    else if (tk == '(')
    {
      //std::cout << "Push OP\n";
      opstack.push(tk);
    }
    else if (tk == ')')
    {
      char top = opstack.top();
      opstack.pop();
      while (top != '(')
      {
        //std::cout << "top: " << top << '\n';
        std::string tops;
        tops += top;
        postfix.push_back(tops);
        top = opstack.top();
        opstack.pop();
      }
    }
    else
    {
      while (! opstack.empty() &&
             (op_prec_char.at(opstack.top()) >= op_prec_char.at(tk)))
      {
        char top = opstack.top();
        std::string tops;
        tops += top;
        postfix.push_back(tops);
        opstack.pop();
      }
      opstack.push(tk);
    }
  }
  while (! opstack.empty())
  {
    char top = opstack.top();
    std::string tops;
    tops += top;
    postfix.push_back(tops);
    opstack.pop();
  }
  return postfix;
}

// Execute the postfix filter stack
double eval_postfix(const std::vector<std::string> pf, const SeidrFileEdge& e)
{
  std::stack<double> stack;
  for (auto s : pf)
  {
    if (op_prec_char.find(s[0]) == op_prec_char.end())
    {
      if (s[0] == 's')
      {
        s.erase(0, 1);
        uint16_t index = std::stoul(s);
        stack.push(e.scores[index].s);
      }
      else if (s[0] == 'r')
      {
        s.erase(0, 1);
        uint16_t index = std::stoul(s);
        stack.push(e.scores[index].r);
      }
      else
      {
        stack.push(std::stod(s));
      }
    }
    else
    {
      double a = stack.top(); stack.pop();
      double b = stack.top(); stack.pop();
      char op = s[0];
      switch (op)
      {
      case '<':
        stack.push(b < a); break;
      case '>':
        stack.push(b > a); break;
      case '=':
        stack.push(almost_equal(a, b)); break;
      case '+':
        stack.push(a + b); break;
      case '-':
        stack.push(a - b); break;
      case '/':
        stack.push(a / b); break;
      case '*':
        stack.push(a * b); break;
      case '|':
        stack.push(a || b); break;
      case '&':
        stack.push(a && b); break;
      }
    }
  }
  double res = stack.top();
  return res;
}

void print_centrality(const seidr_param_view_t& param,
                      const SeidrFileHeader& h,
                      std::ostream& out)
{
  if (param.titles)
  {
    std::vector<std::string> headers = {"Node"};
    if (h.attr.pagerank_calc)
      headers.push_back("PageRank");
    if (h.attr.closeness_calc)
      headers.push_back("Closeness");
    if (h.attr.betweenness_calc)
      headers.push_back("Betweenness");
    if (h.attr.strength_calc)
      headers.push_back("Strength");
    if (h.attr.eigenvector_calc)
      headers.push_back("Eigenvector");
    if (h.attr.katz_calc)
      headers.push_back("Katz");
    for (uint32_t i = 0; i < headers.size(); i++)
    {
      out << headers[i]
          << (i == headers.size() - 1 ? '\n' : '\t');
    }
  }
  h.print_centrality(out);
}

void finish_binary(const seidr_param_view_t& param,
                   SeidrFileHeader& hh)
{
  SeidrFile tf(param.tempfile.c_str());
  tf.open("r");
  SeidrFileHeader th;
  th.unserialize(tf);

  SeidrFile of2(param.ot_file.c_str());
  of2.open("w");
  double ne = hh.attr.edges;
  double nn = hh.attr.nodes;
  if (ne < ((nn * (nn - 1) / 2) * 0.66))
    hh.attr.dense = 0;
  hh.serialize(of2);

  uint32_t i = 1;
  uint32_t j = 0;
  for (uint64_t x = 0; x < hh.attr.edges; x++)
  {
    SeidrFileEdge e;
    e.unserialize(tf, th);

    if (hh.attr.dense)
    {
      while (i != e.index.i && j != e.index.j)
      {
        SeidrFileEdge eo;
        eo.serialize(tf, th);
        j++;
        if (i == j)
        {
          i++;
          j = 0;
        }
      }
      e.serialize(of2, hh);
      j++;
      if (i == j)
      {
        i++;
        j = 0;
      }
    }
    else
    {
      e.serialize(of2, hh);
    }
  }
  of2.close();
  tf.close();
  fs::remove(fs::path(tf.filepath));
}

bool print_edge(const SeidrFileEdge& e, const SeidrFileHeader& h,
                std::ostream& out,
                const std::vector<std::string>& pf,
                const seidr_param_view_t& param,
                SeidrFile& bin_sf, SeidrFileHeader& bin_sf_h,
                bool& do_filter)
{
  bool was_printed = false;
  if (do_filter)
  {
    bool y = ! almost_equal(eval_postfix(pf, e), 0);
    if (y && (! param.binout))
    {
      e.print(out, h, '\n', param.print_supp, param.delim, param.sc_delim,
              param.full, param.tags, param.no_names);
      was_printed = true;
    }
    else if (y && param.binout)
    {
      e.serialize(bin_sf, bin_sf_h);
      bin_sf_h.attr.edges++;
      was_printed = true;
    }
  }
  else
  {
    double x = param.trank ? e.scores[param.tpos].r : e.scores[param.tpos].s;
    bool cut = param.trank ? x < param.threshold : x > param.threshold;
    if (cut && (! param.binout))
    {
      e.print(out, h, '\n', param.print_supp, param.delim, param.sc_delim,
              param.full, param.tags, param.no_names);
      was_printed = true;
    }
    else if (cut && param.binout)
    {
      e.serialize(bin_sf, bin_sf_h);
      bin_sf_h.attr.edges++;
      was_printed = true;
    }
  }
  return was_printed;
}

bool max_lines_reached(const po::variables_map& vm,
                       const uint64_t& nlines,
                       const seidr_param_view_t& param)
{
  if (vm.count("lines"))
  {
    if (nlines >= param.nlines)
    {
      return true;
    }
  }
  return false;
}

bool print_and_check_maxl(const SeidrFileEdge& e, const SeidrFileHeader& h,
                          std::ostream& out,
                          const std::vector<std::string>& pf,
                          const seidr_param_view_t& param,
                          SeidrFile& bin_sf, SeidrFileHeader& bin_sf_h,
                          bool& do_filter,
                          uint64_t& lines_printed,
                          SeidrFile& rf,
                          const po::variables_map& vm)
{
  if (print_edge(e, h, out, pf, param, bin_sf, bin_sf_h, do_filter))
  {
    lines_printed++;
  }
  if (max_lines_reached(vm, lines_printed, param))
  {
    rf.close();
    if (param.binout)
    {
      bin_sf.close();
      finish_binary(param, bin_sf_h);
    }
    return true;
  }
  return false;
}

void print_titles(const seidr_param_view_t& param,
                  const SeidrFileHeader& h,
                  std::ostream& out)
{
  out << "Source\tTarget\tType\t";
  for (uint16_t i = 0; i < h.attr.nalgs; i++)
  {
    out << h.algs[i] << "_score" << param.delim << h.algs[i] << "_rank";
    if (i == h.attr.nalgs - 1 && h.attr.nsupp == 0 && (! param.full))
      out << '\n';
    else
      out << '\t';
  }
  for (uint16_t i = 0; i < h.attr.nsupp; i++)
  {
    out << h.supp[i];
    if (i == h.attr.nsupp - 1 && (! param.full))
      out << '\n';
    else if (i == h.attr.nsupp - 1 && param.full)
      out << '\t';
    else
      out << param.sc_delim;
  }
  if (param.full)
  {
    std::stringstream ost;
    if (h.attr.pagerank_calc)
      ost << "PageRank_Source\tPageRank_target\t";
    if (h.attr.closeness_calc)
      ost << "Closeness_Source\tCloseness_target\t";
    if (h.attr.betweenness_calc)
      ost << "Betweenness_Source\tBetweenness_target\t";
    if (h.attr.strength_calc)
      ost << "Strength_Source\tStrength_target\t";
    if (h.attr.eigenvector_calc)
      ost << "Eigenvector_Source\tEigenvector_target\t";
    if (h.attr.katz_calc)
      ost << "Katz_Source\tKatz_target\t";
    std::string ox = ost.str();
    if (ox.size() > 0)
    {
      ox.back() = '\n';
    } 
    else
    {
      ox += '\n';
    }
    out << ox;
  }
}

int view(int argc, char * argv[]) {

  logger log(std::cerr, "view");

  seidr_param_view_t param;
  uint64_t lines_printed = 0;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " view";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  po::options_description
  umbrella("View and filter contents of SeidrFiles");

  po::options_description opt("Common Options");
  opt.add_options()
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("help,h", "Show this help message")
  ("tempdir,T",
   po::value<std::string>(&param.tempdir)->default_value("", "auto"),
   "Directory to store temporary data");

  po::options_description fopt("Filter Options");
  fopt.add_options()
  ("threshold,t",
   po::value<double>(&param.threshold)->
   default_value(-std::numeric_limits<double>::infinity()),
   "Only print edges with a weight >= t")
  ("threshold-rank,r",
   po::bool_switch(&param.trank)->default_value(false),
   "Threshold edges with a rank of <= t instead")
  ("index,i",
   po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
   "Score column to use as edge weights")
  ("filter,F",
   po::value<std::string>(&param.filter),
   "Supply a filter function to select edges")
  ("nodelist,n",
   po::value<std::string>(&param.nodes),
   "Include only these nodes")
  ("lines,l",
   po::value<uint32_t>(&param.nlines),
   "View only first l results");

  po::options_description Fopt("Formatting Options");
  Fopt.add_options()
  ("no-names,N",
   po::bool_switch(&param.no_names)->default_value(false),
   "Print node index instead of name")
  ("column-headers,c",
   po::bool_switch(&param.titles)->default_value(false),
   "Print column headers")
  ("tags,a",
   po::bool_switch(&param.tags)->default_value(false),
   "Print supplementary information tags")
  ("dense,D",
   po::bool_switch(&param.full)->default_value(false),
   "Print as much information as possible for each edge")
  ("sc-delim,s",
   po::value<std::string>(&param.sc_delim)->default_value(";"),
   "Delimiter for supplementary tags")
  ("delim,d",
   po::value<std::string>(&param.delim)->default_value(";"),
   "Delimiter for scores/ranks");

  po::options_description oopt("Output Options");
  oopt.add_options()
  ("outfile,o",
   po::value<std::string>(&param.ot_file)->default_value("-"),
   "Output file name ['-' for stdout]")
  ("binary,b",
   po::bool_switch(&param.binout)->default_value(false),
   "Output binary SeidrFile");

  po::options_description ocopt("View metadata");
  ocopt.add_options()
  ("header,H",
   po::bool_switch(&param.header)->default_value(false),
   "Print file header as text")
  ("centrality,C",
   po::bool_switch(&param.centrality)->default_value(false),
   "Print node centrality scores if present");

  po::options_description req("Required [can be positional]");
  req.add_options()
  ("in-file", po::value<std::string>(&param.el_file)->required(),
   "Input SeidrFile");

  umbrella.add(req).add(fopt).add(Fopt).add(ocopt).add(oopt).add(opt);

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

  try // Argument sanity
  {
    po::notify(vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  try // General sanity checks
  {
    param.el_file = to_absolute(param.el_file);
    assert_exists(param.el_file);
    assert_can_read(param.el_file);

    if (param.ot_file != "-")
    {
      param.ot_file = to_absolute(param.ot_file);
      if (! param.force)
      {
        assert_no_overwrite(param.ot_file);
      }
      assert_dir_is_writeable(dirname(param.ot_file));
    }

    param.indexfile = param.el_file + ".sfi";
    if (vm.count("nodelist") && (! file_exists(param.indexfile)))
    {
      throw std::runtime_error("SeidrFile index " + param.indexfile +
                               " must exist when using --nodelist");
    }

    assert_mutually_exclusive(vm, {"binary", "header", "centrality", "dense"});
    assert_mutually_exclusive(vm, {"binary", "column-headers"});


    if (param.binout && param.ot_file == "-")
    {
      errno = EINVAL;
      throw std::runtime_error("Writing binary output to stdout is currently "
                               "not supported");
    }
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }


  SeidrFile rf(param.el_file.c_str());
  rf.open("r");
  SeidrFileHeader h;
  h.unserialize(rf);

  std::vector<std::string> pf;
  bool do_filter = false;
  if (vm.count("filter"))
  {
    pf = infix_to_postfix(param.filter);
    do_filter = true;
  }

  make_tpos(param.tpos, h);

  // Prepare binary output structures
  SeidrFileHeader hh = h;
  hh.attr.edges = 0;
  hh.attr.dense = 0;
  param.tempfile = tempfile(param.tempdir);
  SeidrFile of(param.tempfile.c_str());
  if (param.binout)
  {
    of.open("w");
    hh.serialize(of);
  }
  // Prepare textual output structures
  std::shared_ptr<std::ostream> out(&std::cout, [](...) {});
  if (param.ot_file != "-" && (! param.binout))
  {
    out.reset(new std::ofstream(param.ot_file.c_str()));
  }

  // Print only header
  if (param.header)
  {
    h.print(*out);
    return 0;
  }
  // Print centrality measures
  if (param.centrality)
  {
    print_centrality(param, h, *out);
    return 0;
  }
  // Print titles if requested
  if (param.titles)
  {
    print_titles(param, h, *out);
  }
  if (param.nodes != "")
  {
    if (file_exists(param.nodes))
    {
      param.nodelist = read_genes(param.nodes, '\n', '\t');
    }
    else
    {
      param.nodelist = tokenize_delim(param.nodes, ",");
    }

    SeidrFile sfi(param.indexfile.c_str());
    sfi.open("r");
    SeidrFileIndex index;
    index.unserialize(sfi, h);
    sfi.close();
    std::set<offset_t> nodeset = index.get_offset_nodelist(param.nodelist);
    for (const offset_t& o : nodeset)
    {
      rf.seek(o.o);
      SeidrFileEdge e;
      e.unserialize(rf, h);
      e.index.i = o.i;
      e.index.j = o.j;
      if (print_and_check_maxl(e, h, *out, pf, param,
                               of, hh, do_filter, lines_printed, rf, vm))
      {
        return 0;
      }
    }
  }
  else if (h.attr.dense)
  {
    for (uint64_t i = 1; i < h.attr.nodes; i++)
    {
      for (uint64_t j = 0; j < i; j++)
      {
        SeidrFileEdge e;
        e.unserialize(rf, h);
        e.index.i = i;
        e.index.j = j;
        if (EDGE_EXISTS(e.attr.flag))
        {
          if (print_and_check_maxl(e, h, *out, pf, param,
                                   of, hh, do_filter, lines_printed, rf, vm))
          {
            return 0;
          }
        }
      }
    }
  }
  else
  {
    for (uint64_t i = 0; i < h.attr.edges; i++)
    {
      SeidrFileEdge e;
      e.unserialize(rf, h);
      if (print_and_check_maxl(e, h, *out, pf, param,
                               of, hh, do_filter, lines_printed, rf, vm))
      {
        return 0;
      }
    }
  }

  rf.close();

  if (param.binout)
    of.close();

  if (param.binout)
    finish_binary(param, hh);

  return 0;
}
