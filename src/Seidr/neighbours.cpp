// Siedr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <neighbours.h>
#include <BSlogger.h>
// External
#include <vector>
#include <string>
#include <algorithm>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <set>
#include <boost/numeric/conversion/cast.hpp>

int neighbours(int argc, char * argv[])
{

  logger log(std::cerr, "neighbours");

  std::string infile;
  uint16_t n;
  uint32_t tpos;
  std::string nodelist;
  bool trank;
  bool unique;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " neighbours";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Get N neighbours for each node or a list of nodes "
                       "in SeidrFile", ' ', version);

    TCLAP::ValueArg<uint32_t>
    arg_tpos("i", "index", "Index of score to use", false,
             0, "int");
    cmd.add(arg_tpos);

    TCLAP::ValueArg<uint16_t>
    arg_n("n", "neighbours", "Number of neighbours", false,
          10, "int");
    cmd.add(arg_n);

    TCLAP::ValueArg<std::string>
    arg_list("g", "genes", "Comma separated list of nodes of interest", false,
             "", "string");
    cmd.add(arg_list);

    TCLAP::SwitchArg
    arg_r("r", "rank", "Use rank instead of score to sort edges", cmd, false);

    TCLAP::SwitchArg
    arg_u("u", "unique", "Only print unique edges", cmd, false);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true, "", "string");
    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    n = arg_n.getValue();
    tpos = arg_tpos.getValue();
    nodelist = arg_list.getValue();
    trank = arg_r.getValue();
    unique = arg_u.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  std::string indexfile = infile + ".sfi";

  try
  {
    infile = to_canonical(infile);
    indexfile = to_canonical(indexfile);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    log(LOG_ERR) << "Make sure that both the SeidrFile and the file index exist. "
                 << "Run `seidr index` to create an index file\n";
    return errno;
  }

  SeidrFile f(infile.c_str());
  f.open("r");
  SeidrFileHeader h;
  h.unserialize(f);

  if (tpos == 0)
    tpos = h.attr.nalgs - 1;
  else
    tpos--;

  log(LOG_INFO) << "Reading indexfile...\n";
  SeidrFile sfi(indexfile.c_str());
  sfi.open("r");
  SeidrFileIndex index;
  index.unserialize(sfi, h);
  sfi.close();

  if (nodelist != "")
  {
    std::vector<std::string> genes = tokenize_delim(nodelist, ",");
    std::vector<SeidrFileEdge> edges_nu;
    std::set<SeidrFileEdge, sfe_score_sort_class> edges_u;
    for (auto& gene : genes)
    {
      std::vector<offset_t> offsets = index.get_offset_node(gene);
      std::vector<SeidrFileEdge> edges;
      for (auto& offset : offsets)
      {
        f.seek(offset.o);
        SeidrFileEdge e;
        e.unserialize(f, h);
        e.index.i = offset.i;
        e.index.j = offset.j;
        edges.push_back(e);
      }
      if (trank)
        SORT(edges.begin(), edges.end(), sfe_score_sort);
      else
        SORT(edges.begin(), edges.end(), sfe_rank_sort);
      for (uint16_t i = 0; i < n; i++)
      {
        if (unique)
          edges_u.insert(edges[i]);
        else
          edges_nu.push_back(edges[i]);
      }
    }
    if (unique)
    {
      for (auto& e : edges_u)
      {
        e.print(h);
      }
    }
    else
    {
      for (auto& e : edges_nu)
      {
        e.print(h);
      }
    }
  }
  else
  {
    log(LOG_INFO) << "Creating score matrix...\n";
    arma::mat gm = read_network_arma(h, f, -std::numeric_limits<double>::infinity(),
                                     tpos, trank);
    std::vector<offset_t> offsets_nu;
    std::set<offset_t> offsets_u;
    for (arma::uword i = 0; i < gm.n_rows; i++)
    {
      arma::vec v = gm.col(i);

      arma::uvec sort_index;
      if (trank)
        sort_index = arma::sort_index(v);
      else
        sort_index = arma::sort_index(v, "descend");

      for (uint16_t j = 0; j < n; j++)
      {
        uint32_t xi = i;
        uint32_t xj = sort_index[j];
        offset_t off = index.get_offset_pair(h.nodes[xi], h.nodes[xj]);
        if (unique)
          offsets_u.insert(off);
        else
          offsets_nu.push_back(off);
      }
    }
    if (unique)
    {
      for (auto& offset : offsets_u)
      {
        f.seek(offset.o);
        SeidrFileEdge e;
        e.unserialize(f, h);
        e.index.i = offset.i;
        e.index.j = offset.j;
        e.print(h);
      }
    }
    else
    {
      for (auto& offset : offsets_nu)
      {
        f.seek(offset.o);
        SeidrFileEdge e;
        e.unserialize(f, h);
        e.index.i = offset.i;
        e.index.j = offset.j;
        e.print(h);
      }
    }
  }
  return 0;
}