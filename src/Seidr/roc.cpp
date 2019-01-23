/**
 * @file
 * @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
 *
 * @param infs Input files to be aggregated in binary format from el2bin
 * @param genes_file File containing names of genes present in the network
 * @param method Aggregation method: 'borda' implements borda counting
 *               'last' ranks missing edges last, 'ignore' ignores missing.
 * @return int 0 if the function succeeded, an error code otherwise.
 *
 * @section DESCRIPTION
 *
 * This function aggregates any number of networks (represented as binary
 * files by the el2bin function) into a meta-network of scores in [0, 1].
 * If the chosen method is borda counting, we calculate the mean rank of
 * each edge, ranking missing vertices last. The top1 algorithm presents
 * the final score of an edge as the best score of any algorithm. The
 * neural algorithm learns how to best predict new edges from a set of
 * gold standard gene interactions.
 */

// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <roc.h>
#include <BSlogger.h>
// External
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include <map>
#include <boost/numeric/conversion/cast.hpp>

using boost::numeric_cast;

std::vector< std::pair< uint32_t, uint32_t > >
make_gold_vec(std::string infile, std::map<std::string, uint32_t>& refmap)
{
  logger log(std::cerr, "ROC::make_gold_vec");
  std::vector< std::pair< uint32_t, uint32_t > > ret;
  std::ifstream ifs(infile.c_str(), std::ios::in);
  uint64_t not_found = 0;
  std::string line;
  while (std::getline(ifs, line))
  {
    std::stringstream ss(line);
    std::string x, y;
    ss >> x; ss >> y;

    auto ptr_a = refmap.find(x);
    auto ptr_b = refmap.find(y);

    if (ptr_a == refmap.end())
    {
      not_found++;
      continue;
    }
    else if (ptr_b == refmap.end())
    {
      not_found++;
      continue;
    }

    uint32_t a = refmap[x];
    uint32_t b = refmap[y];

    if (a < b)
    {
      uint32_t c = a;
      a = b;
      b = c;
    }
    ret.push_back(std::pair<uint32_t, uint32_t>(a, b));
  }
  if (not_found > 0)
  {
    log(LOG_WARN) << not_found << " gold standard edges had at least one "
                  << "node which was not present in the network and were "
                  << "omitted.\n";
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

std::vector<MiniEdge> filter_tf(std::string tfs,
                                std::vector<MiniEdge>& nv,
                                std::map<std::string, uint32_t>& refmap)
{
  std::vector<MiniEdge> ret;
  std::vector<uint32_t> tfs_mapped;
  std::ifstream ifs(tfs.c_str());

  std::string line;
  while (std::getline(ifs, line))
    tfs_mapped.push_back(refmap[line]);

  std::sort(tfs_mapped.begin(), tfs_mapped.end());
  for (auto& e : nv)
  {
    if (std::binary_search(tfs_mapped.begin(), tfs_mapped.end(), e.i) ||
        std::binary_search(tfs_mapped.begin(), tfs_mapped.end(), e.j))
      ret.push_back(e);
  }

  return ret;
}

std::vector<double> make_range_cutoffs(std::vector< MiniEdge >& v, uint32_t n)
{
  double x = n;
  double min = v[v.size() - 1].s;
  double max = v[0].s;
  double alpha = (max - min) / x;
  std::vector< double > ret;
  for (uint32_t i = 0; i < n; i++)
  {
    max -= alpha;
    ret.push_back(max);
  }
  return ret;
}

std::pair< uint32_t, uint32_t >
get_tp_fp(std::vector< std::pair< uint32_t, uint32_t> >& truth,
          std::vector< MiniEdge >& v)
{
  uint32_t tp = 0, fp = 0;
  for (auto& e : v)
  {
    uint32_t i = e.i;
    uint32_t j = e.j;
    if (i < j)
    {
      uint32_t k = i;
      i = j;
      j = k;
    }
    auto x = std::pair<uint32_t, uint32_t>(i, j);
    if (std::binary_search(truth.begin(), truth.end(), x))
      tp++;
    else
      fp++;
  }
  return std::pair<uint32_t, uint32_t>(tp, fp);
}

void resize_by_fraction(std::vector< std::pair< uint32_t, uint32_t> >& truth,
                        std::vector< MiniEdge >& v,
                        double& fedg)
{
  logger log(std::cerr, "ROC::resize_by_fraction");
  uint32_t tp = 0, fp = 0;
  uint64_t cnt = 0;
  for (auto& e : v)
  {
    uint32_t i = e.i;
    uint32_t j = e.j;
    if (i < j)
    {
      uint32_t k = i;
      i = j;
      j = k;
    }
    auto x = std::pair<uint32_t, uint32_t>(i, j);
    if (std::binary_search(truth.begin(), truth.end(), x))
      tp++;
    else
      fp++;

    cnt++;

    double f = numeric_cast<double> (tp) / numeric_cast<double> (truth.size());
    if (f >= fedg)
    {
      log(LOG_INFO) << "Resizing to top " << cnt << " egdes in order to keep "
      << "at least " << truth.size() * f << " TP edges.\n";
      v.resize(cnt);
    }
  }

}

void print_roc(std::vector<std::pair<uint32_t, uint32_t>>& truth,
               std::vector<MiniEdge>& v,
               std::pair<uint32_t, uint32_t>& tpfp,
               uint32_t& datap)
{
  logger log(std::cerr, "ROC");
  uint32_t cind = 0;
  uint32_t tp = 0, fp = 0;
  uint32_t cnt = 0;
  uint32_t ne = v.size();
  uint32_t intervx = 0;
  uint32_t intervy = 0;
  if (datap != 0)
  {
    intervy = numeric_cast<double> (ne) /
              numeric_cast<double> (datap);
  }
  log(LOG_INFO) << "intervy: " << intervy << '\n';
  std::cout << "#TP/FP\t" <<  tpfp.first << '\t' << tpfp.second << '\n';
  for (auto& e : v)
  {
    cnt++;
    uint32_t i = e.i;
    uint32_t j = e.j;
    if (i < j)
    {
      uint32_t k = i;
      i = j;
      j = k;
    }
    auto x = std::pair<uint32_t, uint32_t>(i, j);
    if (std::binary_search(truth.begin(), truth.end(), x))
      tp++;
    else
      fp++;
    double tpr = numeric_cast<double> (tp) / numeric_cast<double> (tpfp.first);
    double fpr = numeric_cast<double> (fp) / numeric_cast<double> (tpfp.second);
    double ppv = numeric_cast<double> (tp) / numeric_cast<double> (cnt);
    if (cnt > intervx || cnt == ne - 1)
    {
      std::cout << tpr << '\t' << fpr << '\t' << ppv << '\n';
      intervx += intervy;
    }
  }
}

int roc(int argc, char * argv[])
{

  logger log(std::cerr, "ROC");

  std::string gold;
  std::string netw;
  uint32_t tpos;
  uint32_t nedg;
  std::string tfs;
  uint32_t datap;
  double fedg;

  const char * args[argc - 1];
  args[0] = strcat(argv[0], " roc");
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Calculate ROC curves of predictions in "
                       "SeidrFiles given true edges", ' ', version);

    TCLAP::ValueArg<std::string>
    arg_gold("g", "gold", "Gold standard (true edges) input file", true,
             "", "string");
    cmd.add(arg_gold);

    TCLAP::ValueArg<std::string>
    arg_network("n", "network", "Network input file", true,
                "", "string");
    cmd.add(arg_network);

    TCLAP::ValueArg<uint32_t>
    arg_tpos("i", "index", "Index of score to use (1-based)", false,
             0, "int");
    cmd.add(arg_tpos);

    TCLAP::ValueArg<uint32_t>
    arg_nedg("e", "edges", "Number of top edges to consider", false,
             0, "int");
    cmd.add(arg_nedg);

    TCLAP::ValueArg<double>
    arg_fedg("E", "fraction", "Fraction of gold standard edges to include", false,
             -1, "int");
    cmd.add(arg_fedg);

    TCLAP::ValueArg<uint32_t>
    arg_datap("p", "points", "Number of data points to print", false,
              0, "int");
    cmd.add(arg_datap);

    TCLAP::ValueArg<std::string>
    arg_tfs("t", "tfs", "List of transcription factors to consider", false,
            "", "string");
    cmd.add(arg_tfs);

    // Parse arguments
    cmd.parse(argc, args);
    gold = arg_gold.getValue();
    netw = arg_network.getValue();
    tpos = arg_tpos.getValue();
    nedg = arg_nedg.getValue();
    tfs = arg_tfs.getValue();
    datap = arg_datap.getValue();
    fedg = arg_fedg.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }
  try
  {
    if (! file_can_read(gold.c_str()))
      throw std::runtime_error("Cannot read file: " + gold);
    if (! file_can_read(netw.c_str()))
      throw std::runtime_error("Cannot read file: " + netw);
    if (tfs != "")
      if (! file_can_read(tfs.c_str()))
        throw std::runtime_error("Cannot read file: " + tfs);

    if (nedg > 0 && fedg > 0)
      throw std::runtime_error("Please provide only one of -e and -E");

    if (fedg > 1)
      throw std::runtime_error("-E should be in (0, 1)");


    SeidrFile f(netw.c_str());
    f.open("r");
    SeidrFileHeader h;
    h.unserialize(f);

    if (tpos > h.attr.nalgs)
    {
      throw std::runtime_error("Index (-i) selected is larger than the number"
                               " of algorithms in the network.");
    }
    if (tpos == 0)
      tpos = h.attr.nalgs - 1;
    else
      tpos--;

    std::map<std::string, uint32_t> refmap;
    uint32_t ind = 0;
    for (auto& no : h.nodes)
    {
      refmap[no] = ind++;
    }

    auto gv = make_gold_vec(gold, refmap);
    auto nv = read_network_minimal(h, f, -1, tpos);

    std::sort(nv.begin(), nv.end(), me_score_sort_abs); //TODO: sort on rank, this is not robust
    if (tfs != "")
      nv = filter_tf(tfs, nv, refmap);
    if (nedg != 0)
      nv.resize(nedg);

    if (fedg > 0)
      resize_by_fraction(gv, nv, fedg);

    auto tpfp = get_tp_fp(gv, nv);

    print_roc(gv, nv, tpfp, datap);

    f.close();
    return 0;
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

}
