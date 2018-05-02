/**
* @file
* @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
* @version 0.01
*
* @section DESCRIPTION
*
* This function takes an edge list of scores and transforms it
* into a binary file of ranks (using sports ranking, e.g: if
* two ranks were tied for the first place, both would assume
* the first rank, but there would be no second rank).
*/

// Seidr
#include <common.h>
#include <import.h>
#include <gzstream.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
//#include <functional>
#include <map>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <memory>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

using boost::lexical_cast;
using boost::numeric_cast;
typedef std::map<std::string, uint32_t> stringmap;
typedef std::pair<uint32_t, uint32_t> inx_t;

read_logger& read_logger::operator++()
{
  _i++;
  if (_i % _n == 0)
  {
    _l(LOG_INFO) << "Read " << _i << " edges\n";
  }
  return *this;
}

read_logger read_logger::operator++(int)
{
  read_logger tmp = *this;
  ++*this;
  return tmp;
}

void lt_map::insert(inx_t& inx, reduced_edge& ee, bool& rev,
                    bool& drop)
{
  // Always work in lower triangle
  if (inx.first < inx.second)
  {
    uint32_t tmp = inx.first;
    inx.first = inx.second;
    inx.second = tmp;
    // Keep track of edge direction
    ee.d = 2;
  }
  else if (ee.d != 3)
  {
    ee.d = 1;
  }
  else
  {
    ee.d = 0;
  }
  edge e;
  e.i = inx.first;
  e.j = inx.second;
  e.d = ee.d;
  e.w = ee.w;
  // Check if edge already exists
  auto it = _ev.find(e);
  // If not, insert edge and do nothing else
  if (it == _ev.end())
  {
    if (drop)
    {
      if (! almost_equal(e.w, 0))
        _ev.insert(e);
    }
    else
    {
      _ev.insert(e);
    }
  }
  // If there is a symmetric edge present, check if current edge has a better
  // score and insert if it has.
  else if ( !drop || (drop && ! almost_equal(e.w, 0)))
  {
    // In reverse mode, bigger is better
    if (rev)
    {
      if (e.w > it->w)
      {
        (*it).w = e.w;
        (*it).d = e.d;
      }
      else if (almost_equal(e.w, it->w))
      {
        // If A->B == B<-A, set to undirected.
        (*it).d = 0;
      }
    }
    // In default mode, lower is better
    else
    {
      if (e.w < it->w)
      {
        (*it).w = e.w;
        (*it).d = e.d;
      }
      else if (almost_equal(e.w, it->w))
      {
        // If A->B == B<-A, set to undirected.
        (*it).d = 0;
      }
    }
  }
}

// Convert lower triangular map to vector of edges
std::vector<const edge*> lt_map::to_vec()
{
  std::vector<const edge*> ret;
  for (auto it = _ev.begin(); it != _ev.end(); it++)
  {
    ret.push_back(&(*it));
    //_ev.erase(it);
  }
  return ret;
}

/**
* A comparison to rank scores in ascending or descending order
*/
bool ascending(const edge* a, const edge* b) { return (a->w < b->w);}
bool descending(const edge* a, const edge* b) { return (a->w > b->w);}
bool abs_ascending(const edge* a, const edge* b) { return (fabs(a->w) < fabs(b->w));}
bool abs_descending(const edge* a, const edge* b) { return (fabs(a->w) > fabs(b->w));}
bool lt(const edge* a, const edge* b)
{
  if (a->i < b->i)
    return true;
  if (a->i == b->i)
    return a->j < b->j;
  return false;
}
/**
* A function to replace the scores in a vector with the
* ranks they assume in the entire set.
*
* @param a A vector iterator pointing to the start of the vector.
* @param b A vector iterator pointing to the end of the vector.
* @param reverse Boolean indicator to check if the vector should
*                be sorted in descending order.
*/
void rank_vector(std::vector<const edge*>& ev, bool reverse, bool absolute)
{
//Sorting
  logger log(std::cerr, "rank_vector");
  if (reverse)
  {
    if (absolute)
      SORT(ev.begin(), ev.end(), abs_descending);
    else
      SORT(ev.begin(), ev.end(), descending);
  }
  else
  {
    if (absolute)
      SORT(ev.begin(), ev.end(), abs_ascending);
    else
      SORT(ev.begin(), ev.end(), ascending);
  }
  auto it = ev.begin();
  uint64_t pos = 0;
  seidr_score_t prev = (*it)->w;
  uint64_t start = 0;
  seidr_score_t rank;
  while (it != ev.end())
  {
    it++; pos++;
    if (it == ev.end() || (*it)->w != prev)
    {
      rank = ( lexical_cast<seidr_score_t>(pos) + 1 +
               lexical_cast<seidr_score_t>(start) ) / 2;
      for (uint64_t i = start; i < pos; i++)
      {
        //edge e = ev[i];
        //e.r = rank;
        ev[i]->r = rank;
      }
      if (it != ev.end())
      {
        start = pos;
        prev = (*it)->w;
      }
    }
  }
  SORT(ev.begin(), ev.end(), lt);
}


/**
* Creating a binary ranked file from an edge list input.
*
* @param el_file The input edge list.
* @param gene_file The gene names.
* @return int 0 if the function succeeded, an error code otherwise.
*/
int import(int argc, char * argv[]) {

  logger log(std::cerr, "import");
  read_logger pr(log); // processed rows

// Variables used by the function
  std::string el_file;
  std::string gene_file;
  std::string out_file;
  bool rev;
  std::vector<std::string> genes;
  std::string name;
  bool is_lm = false;
  bool is_aracne = false;
  bool is_full_matrix = false;
  bool is_edge_list = false;
  bool force_undirected;
  bool absolute;
  bool drop_zero;
  bool force = false;
  std::string format;

// We ignore the first argument
  const char * args[argc - 1];
  std::string p(argv[0]);
  p += " import";
  args[0] = p.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;


  try
  {
// Add arguments from the command line
    TCLAP::CmdLine cmd("Convert various text base network representations to"
                       " SeidrFiles", ' ', version);

    TCLAP::ValueArg<std::string>
    arg_genes("g", "genes", "File containing gene names",
              true, "", "");
    cmd.add(arg_genes);

    std::vector<std::string> type_constraints{"el", "lm", "m", "ara"};
    TCLAP::ValuesConstraint<std::string> t_constraints(type_constraints);
    TCLAP::ValueArg<std::string>
    arg_format("f", "format", "Input file format.", false,
               "el", &t_constraints);
    cmd.add(arg_format);

    TCLAP::ValueArg<std::string>
    arg_infile("i", "infile", "Input file. '-' for stdin", true,
               "", "");
    cmd.add(arg_infile);

    TCLAP::ValueArg<std::string>
    arg_outfile("o", "outfile", "Output file", false,
                "elranks.sf", "");
    cmd.add(arg_outfile);

    TCLAP::ValueArg<std::string>
    arg_name("n", "name", "Name of the network", false,
             "unnamed", "unnamed");

    cmd.add(arg_name);

// Switches
    TCLAP::SwitchArg 
    arg_absolute("A", "absolute", "Rank absolutes of values.", cmd, false);

    TCLAP::SwitchArg 
    arg_reverse("r", "reverse", "Higher scores/values are better.", cmd, false);

    TCLAP::SwitchArg 
    arg_drop_zero("z", "drop_zero", 
                  "Drop edges with scores of zero from the network", cmd, false);

    TCLAP::SwitchArg 
    arg_undirected("u", "undirected", 
                   "Force edges to be interpreted as undirected.", cmd, false);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

// Parse arguments
    cmd.parse(argc, args); // FIX in deploy
    gene_file = arg_genes.getValue();
    el_file = arg_infile.getValue();
    rev = arg_reverse.getValue();
    out_file = arg_outfile.getValue();
    name = arg_name.getValue();
    absolute = arg_absolute.getValue();
    drop_zero = arg_drop_zero.getValue();
    format = arg_format.getValue();
    force_undirected = arg_undirected.getValue();
    force = switch_force.getValue();

    if (format == "el")
      is_edge_list = true;
    else if (format == "lm")
      is_lm = true;
    else if (format == "m")
      is_full_matrix = true;
    else if (format == "ara")
      is_aracne = true;

  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << "[Invalid argument exception]: " << except.what() << '\n';
    return EINVAL;
  }

  try
  {
    if (el_file != "-")
    {
      if (! file_can_read(el_file.c_str()))
        throw std::runtime_error("Cannot read file " + el_file);
    }
    if (! file_can_read(gene_file.c_str()))
    {
      throw std::runtime_error("Cannot read file " + gene_file);
    }
    if (! file_can_create(out_file.c_str()))
    {
      throw std::runtime_error("Cannot create file " + out_file);
    }
    if (is_lm + is_aracne + is_full_matrix + is_edge_list != 1)
    {
      throw std::runtime_error("Exactly one input format must be specified");
    }
    if (file_exists(out_file) && ! force)
      throw std::runtime_error("File exists: " + out_file);
  }
  catch (std::exception& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return errno;
  }

  genes = read_genes(gene_file, '\n', '\t');
  stringmap g_map;

  size_t i = 0;
  for (std::string& s : genes)
  {
    g_map.insert(std::pair<std::string, size_t>(s, i++));
  }

  log(LOG_INFO) << "Read " << genes.size() << " gene names\n";

  log(LOG_INFO) << "Allocating initial matrix\n";
  lt_map X;

  X._ev.reserve((genes.size() * (genes.size() - 1)) / 2);

  log(LOG_INFO) << "Populating matrix\n";
  std::shared_ptr<std::istream> ifs;

  if (el_file == "-")
  {
    ifs.reset(&(std::cin), [](...) {});
  }
  else
  {
    ifs.reset(new std::ifstream(el_file.c_str(), std::ifstream::in));
  }

  if (is_edge_list)
  {
    try
    {
//    std::istream ifs = *in_stream;
      if (! (*ifs).good())
      {
        throw std::runtime_error("Error opening edge list. Check permissions/spelling.");
      }
      seidr_score_t v = 0; // value
      std::string gi;
      std::string gj;

      for (std::string row; std::getline((*ifs), row, '\n');)
      {
        std::istringstream fl(row);
        size_t x = 0;
        for (std::string field; std::getline(fl, field, '\t');)
        {
          if (x == 0) gi = field;
          if (x == 1) gj = field;
#ifndef SEIDR_SCORE_DOUBLE
          if (x == 2) v = std::stof(field);
#else
          if (x == 2) v = std::stod(field);
#endif
          x++;
        }
        reduced_edge e;
        inx_t inx;
        inx.first = g_map[gi];
        inx.second = g_map[gj];
        e.w = v;
        if (force_undirected)
          e.d = 3;
        X.insert(inx, e, rev, drop_zero);
        pr++;
      }
      if ((*ifs).bad())
      {
        throw std::runtime_error("Error reading edge list.");
      }
//(*ifs).close();
    }
    catch (std::runtime_error& except)
    {
      log(LOG_ERR) << "[Runtime error]: " << except.what() << '\n';
      return errno;
    }
  }
  else if (is_lm)
  {
//Lower triangular matrix, no diagonal
//      std::istream ifs = *in_stream;
    if (! (*ifs).good())
    {
      throw std::runtime_error("Error opening edge list. Check permissions/spelling.");
    }
    for (uint32_t i = 1; i < genes.size(); i++)
    {
      std::string line;
      std::getline((*ifs), line);
      std::stringstream ss(line);
      for (uint32_t j = 0; j < i; j++)
      {
        std::string col;
        ss >> col;
#ifndef SEIDR_SCORE_DOUBLE
        seidr_score_t x = std::stof(col);
#else
        seidr_score_t x = std::stod(col);
#endif
        reduced_edge e;
        inx_t inx;
        inx.first = i;
        inx.second = j;
        e.w = x;
        if (force_undirected)
          e.d = 3; // Force undirected
        X.insert(inx, e, rev, drop_zero);
        if (pr._i % 1000000 == 0)
        {
          log(LOG_INFO) << "BC: " << X._ev.bucket_count()
                        << ", LF: " << X._ev.load_factor() << '\n';
        }
        pr++;
      }
    }
  }
  else if (is_full_matrix)
  {
    if (! (*ifs).good())
    {
      throw std::runtime_error("Error opening edge list. Check permissions/spelling.");
    }
    for (uint32_t i = 0; i < genes.size(); i++)
    {
      std::string line;
      std::getline((*ifs), line);
      std::stringstream ss(line);
      for (uint32_t j = 0; j < genes.size(); j++)
      {
        std::string col;
        ss >> col;
        if (i == j) continue;
#ifndef SEIDR_SCORE_DOUBLE
        seidr_score_t x = std::stof(col);
#else
        seidr_score_t x = std::stod(col);
#endif
        reduced_edge e;
        inx_t inx;
        inx.first = i;
        inx.second = j;
        e.w = x;
        if (force_undirected)
          e.d = 3;
        X.insert(inx, e, rev, drop_zero);
        pr++;
      }
    }
  }
  else if (is_aracne)
  {
//ARACNE code
//      std::istream ifs = *in_stream;
    if (! (*ifs).good())
    {
      throw std::runtime_error("Error opening edge list. Check permissions/spelling.");
    }
    std::string line;
    while (std::getline( (*ifs), line))
    {
      if (line.at(0) == '>') continue; //comments
      std::stringstream ss(line);
      std::string gi;
      ss >> gi;
      uint32_t i = g_map.at(gi);
      uint32_t tmp = 1;
      uint32_t j;
      seidr_score_t x;
      std::string col;
      std::string target;
      while (ss.good())
      {
        ss >> col;
        if (tmp % 2 == 1)
        {
          target = col;
        }
        else if (tmp % 2 == 0)
        {
          j = g_map.at(target);
#ifndef SEIDR_SCORE_DOUBLE
          x = std::stof(col);
#else
          x = std::stod(col);
#endif
          reduced_edge e;
          inx_t inx;
          inx.first = i;
          inx.second = j;
          e.w = x;
          if (force_undirected)
            e.d = 3;
          X.insert(inx, e, rev, drop_zero);
          pr++;
        }
        tmp++;
      }
    }
  }
  else
  {
//Shouldn't happen...
    throw std::runtime_error("Unknown file format.");
  }

  //log(LOG_INFO) << "Read " << pr << " edges\n";

  log(LOG_INFO) << "Swapping min/max edges to lower triangle\n";
  std::vector<const edge*> vs = X.to_vec();

  log(LOG_INFO) << "Computing ranks\n";
  rank_vector(vs, rev, absolute);


  log(LOG_INFO) << "Writing data: "
                << genes.size() << " edges, "
                << vs.size()  << " nodes\n";

  try
  {
    uint64_t cnt = 0;
    uint8_t dense = 0;
    double nn = genes.size();
    double ne = vs.size();
    if (ne > ((nn * (nn - 1) / 2) * 0.66))
      dense = 1;

    SeidrFile ostr(out_file.c_str());
    ostr.open("w");

    SeidrFileHeader h;
    h.attr.nodes = genes.size();
    h.attr.edges = vs.size();
    h.attr.nalgs = 1;
    h.attr.dense = dense;
    h.version_from_char(_XSTR(VERSION));
    h.cmd_from_args(argv, argc);


    for (auto it = genes.begin(); it != genes.end(); it++)
    {
      h.nodes.push_back(*it);
    }

    std::vector<std::string> method = {name};
    h.algs = method;
    h.serialize(ostr);

// Main data. di and dj track the lower triangular index
// of the current edge and empty edges will be created
// in case they are missing in the rank vector. In sparse
// mode these are ignored
    uint32_t di = 1;
    uint32_t dj = 0;
    uint64_t end = dense ? (genes.size() * (genes.size() - 1)) / 2 : vs.size();
    for (size_t i = 0; i < end; i++)
    {
      SeidrFileEdge e;

// Create an empty entry for missing edges in dense
// storage mode
      if (dense && i < vs.size())
      {
        while (vs[i]->i != di || vs[i]->j != dj)
        {
          SeidrFileEdge tmp;
          EDGE_SET_MISSING(tmp.attr.flag);
          tmp.serialize(ostr, h);
          if (dj == di - 1)
          {
            di++;
            dj = 0;
          }
          else
          {
            dj++;
          }
        }
// If i,j of the lower triangular index are euqal
// to the current node, increment them to the next
// node in the index for the next iteration
        if (dj == di - 1)
        {
          di++;
          dj = 0;
        }
        else
        {
          dj++;
        }
      }

// Only relevant in sparse mode
      if (! dense)
      {
        e.index.i = vs[i]->i;
        e.index.j = vs[i]->j;
      }
      if (i < vs.size())
      {
        edge_score s;
        s.s = vs[i]->w;
        s.r = vs[i]->r;
        e.scores.push_back(s);

        EDGE_SET_EXISTING(e.attr.flag);
        switch (vs[i]->d)
        {
        case 0:
          EDGE_SET_UNDIRECTED(e.attr.flag);
          break;
        case 1:
          EDGE_SET_AB(e.attr.flag);
          break;
        case 2:
          EDGE_SET_BA(e.attr.flag);
          break;
        default:
          throw std::runtime_error("Unknown edge direction\n");
        }
        e.serialize(ostr, h);
        cnt++;
        if (cnt % 100000 == 0)
        {
          log(LOG_INFO) << cnt << " edges written ("
                        << (numeric_cast<double>(cnt) / numeric_cast<double>(vs.size())) * 100.0
                        << "%)\n";
        }
      }
      else
      {
        SeidrFileEdge tmp;
        EDGE_SET_MISSING(tmp.attr.flag);
        tmp.serialize(ostr, h);
      }
    }
    ostr.close();
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
  }

  log.time_since_start();
  return 0;
}
