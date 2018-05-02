/**
 * @file
 * @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
 * @version 0.01
 *
 * @param el_file Binary representation of ranked edge list or aggregated network.
 * @param gene_file File with gene names.
 * @return int 0 if the function succeeded, an error code otherwise.
 *
 * @section DESCRIPTION
 *
 * This function simply converts the binary edge lists into human readble text.
 */
// Seidr
#include <common.h>
#include <viewRanks.h>
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
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>
#include <map>

namespace fs = boost::filesystem;

const static std::map<char, uint8_t> op_prec_char =
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
double eval_postfix(std::vector<std::string> pf, SeidrFileEdge& e)
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

void print_edgelist(seidr_param_view_t& param)
{

}

int view(int argc, char * argv[]) {

  logger log(std::cerr, "view");

  seidr_param_view_t param;
  param.print_supp = true;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " view";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("View or filter SeidrFiles.", ' ', version);

    TCLAP::ValueArg<double>
    arg_threshold("t", "threshold", "Only print edges with a weight >= t", false,
                  -std::numeric_limits<double>::infinity(), "-Inf");
    cmd.add(arg_threshold);

    TCLAP::ValueArg<uint32_t>
    arg_tindex("i", "index", "Apply threshold on this index", false,
               0, "float");
    cmd.add(arg_tindex);

    TCLAP::ValueArg<std::string>
    arg_filter("f", "filter", "Supply a filter function to select edges", false,
               "", "string");
    cmd.add(arg_filter);

    TCLAP::ValueArg<std::string>
    arg_nodes("n", "nodelist", "Filter specific nodes", false,
              "", "string");
    cmd.add(arg_nodes);

    TCLAP::SwitchArg
    switch_no_names("N", "no-names", "Do not resolve node names"
                    " (print index instead)", cmd, false);

    TCLAP::ValueArg<std::string>
    arg_outfile("o", "outfile", "Output to file", false,
                "", "string");
    cmd.add(arg_outfile);

    TCLAP::ValueArg<std::string>
    arg_scdelim("s", "sc-delim", "Delimiter for supplementary tags", false,
                ";", "string");
    cmd.add(arg_scdelim);

    TCLAP::ValueArg<std::string>
    arg_delim("d", "delim", "Delimiter for scores/ranks", false,
              ";", "string");
    cmd.add(arg_delim);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file (aggregated gene counts)", true,
               "", "string");
    cmd.add(arg_infile);

    TCLAP::SwitchArg
    arg_h("H", "header", "Print header of file", cmd, false);

    TCLAP::SwitchArg
    arg_t("T", "title", "Print column titles", cmd, false);

    TCLAP::SwitchArg
    arg_c("C", "centrality", "Output centrality scores", cmd, false);

    TCLAP::SwitchArg
    arg_f("F", "full", "Output full information", cmd, false);

    TCLAP::SwitchArg
    arg_r("r", "threshold-rank",
          "Apply threshold value to rank rather than score", cmd, false);

    TCLAP::SwitchArg
    arg_b("b", "binary", "Output binary format", cmd, false);

    TCLAP::SwitchArg
    arg_g("a", "tags", "Print supplementary information tags", cmd, false);
    // Parse arguments
    cmd.parse(argc, args);
    // Get parameters
    param.el_file = arg_infile.getValue();
    param.header = arg_h.getValue();
    param.titles = arg_t.getValue();
    param.threshold = arg_threshold.getValue();
    param.tpos = arg_tindex.getValue();
    param.trank = arg_r.getValue();
    param.filter = arg_filter.getValue();
    param.ot_file = arg_outfile.getValue();
    param.binout = arg_b.getValue();
    param.centrality = arg_c.getValue();
    param.full = arg_f.getValue();
    param.tags = arg_g.getValue();
    param.sc_delim = arg_scdelim.getValue();
    param.delim = arg_delim.getValue();
    param.nodes = arg_nodes.getValue();
    param.no_names = switch_no_names.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    errno = EINVAL;
    log(LOG_ERR) << except.what() << '\n';
    return errno;
  }

  uint16_t incomp_args = 0;
  incomp_args += param.header;
  incomp_args += param.binout;
  incomp_args += param.centrality;
  incomp_args += param.full;

  if (incomp_args > 1)
  {
    errno = EINVAL;
    throw std::runtime_error("Incompatible Arguments: Only supply one of "
                             "--header, --binary, --centrality or --full");
  }

  if (! file_can_read(param.el_file.c_str()))
  {
    errno = ENOENT;
    throw std::runtime_error("Cannot read file: " + param.el_file);
  }

  SeidrFile rf(param.el_file.c_str());
  rf.open("r");

  SeidrFileHeader h;
  h.unserialize(rf);

  std::vector<std::string> pf;
  bool do_filter = false;
  if (param.filter != "")
  {
    pf = infix_to_postfix(param.filter);
    do_filter = true;
  }

  if (param.tpos == 0)
    param.tpos = h.attr.nalgs - 1;
  else
    param.tpos--;

  SeidrFileHeader hh = h;
  hh.attr.edges = 0;
  hh.attr.dense = 0;

  std::shared_ptr<std::ostream> out(&std::cout, [](...) {});

  if (param.ot_file != "" && (! param.binout))
  {
    out.reset(new std::ofstream(param.ot_file.c_str()));
  }
  else if (param.ot_file == "" && param.binout)
  {
    param.ot_file = "thresholded.sf";
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
        (*out) << headers[i]
               << (i == headers.size() - 1 ? '\n' : '\t');
      }
    }
    h.print_centrality(*out);
    return 0;
  }

  std::string tmpfile = param.ot_file + ".tmp";
  SeidrFile of(tmpfile.c_str());
  if (param.binout)
  {
    of.open("w");
    hh.serialize(of);
  }

  if (param.titles && (! param.binout))
  {
    *out << "Source\tTarget\tType\t";
    for (uint16_t i = 0; i < h.attr.nalgs; i++)
    {
      *out << h.algs[i] << "_score" << param.delim << h.algs[i] << "_rank";
      if (i == h.attr.nalgs - 1 && h.attr.nsupp == 0 && (! param.full))
        *out << '\n';
      else
        *out << '\t';
    }
    for (uint16_t i = 0; i < h.attr.nsupp; i++)
    {
      *out << h.supp[i];
      if (i == h.attr.nsupp - 1 && (! param.full))
        *out << '\n';
      else if (i == h.attr.nsupp - 1 && param.full)
        *out << '\t';
      else
        *out << param.sc_delim;
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
      ox.back() = '\n';
      *out << ox;
    }
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

    std::string indexfile = param.el_file + ".sfi";
    if (! file_exists(indexfile))
      throw std::runtime_error("SeidrFile index " + indexfile + " must exist "
                               "when using --nodelist");

    SeidrFile sfi(indexfile.c_str());
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
      if (do_filter)
      {
        bool y = ! almost_equal(eval_postfix(pf, e), 0);
        if (y && (! param.binout))
          e.print(*out, h, '\n', param.print_supp, param.delim, param.sc_delim,
                  param.full, param.tags, param.no_names);
        else if (y && param.binout)
        {
          e.serialize(of, hh);
          hh.attr.edges++;
        }
      }
      else
      {
        double x = param.trank ? e.scores[param.tpos].r : e.scores[param.tpos].s;
        bool cut = param.trank ? x < param.threshold : x > param.threshold;
        if (cut && (! param.binout))
          e.print(*out, h, '\n', param.print_supp, param.delim, param.sc_delim,
                  param.full, param.tags, param.no_names);
        else if (cut && param.binout)
        {
          e.serialize(of, hh);
          hh.attr.edges++;
        }
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
          if (do_filter)
          {
            bool y = ! almost_equal(eval_postfix(pf, e), 0);
            if (y && (! param.binout))
              e.print(*out, h, '\n', param.print_supp, param.delim, param.sc_delim,
                      param.full, param.tags, param.no_names);
            else if (y && param.binout)
            {
              e.serialize(of, hh);
              hh.attr.edges++;
            }
          }
          else
          {
            double x = param.trank ? e.scores[param.tpos].r : e.scores[param.tpos].s;
            bool cut = param.trank ? x < param.threshold : x > param.threshold;
            if (cut && (! param.binout))
              e.print(*out, h, '\n', param.print_supp, param.delim, param.sc_delim,
                      param.full, param.tags, param.no_names);
            else if (cut && param.binout)
            {
              e.serialize(of, hh);
              hh.attr.edges++;
            }
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
      if (do_filter)
      {
        bool y = ! almost_equal(eval_postfix(pf, e), 0);
        if (y && (! param.binout))
          e.print(*out, h, '\n', param.print_supp, param.delim, param.sc_delim,
                  param.full, param.tags, param.no_names);
        else if (y && param.binout)
        {
          e.serialize(of, hh);
          hh.attr.edges++;
        }
      }
      else
      {
        double x = param.trank ? e.scores[param.tpos].r : e.scores[param.tpos].s;
        bool cut = param.trank ? x < param.threshold : x > param.threshold;
        if (cut && (! param.binout))
          e.print(*out, h, '\n', param.print_supp, param.delim, param.sc_delim,
                  param.full, param.tags, param.no_names);
        else if (cut && param.binout)
        {
          e.serialize(of, hh);
          hh.attr.edges++;
        }
      }
    }
  }
  rf.close();

  if (param.binout)
    of.close();

  if (param.binout)
  {
    SeidrFile tf(tmpfile.c_str());
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

  return 0;
}
