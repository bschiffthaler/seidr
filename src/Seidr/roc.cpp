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
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <roc.h>
#include <BSlogger.h>
#include <parallel_control.h>
// Parallel includes
#if defined(SEIDR_PSTL)
#include <tbb/task_scheduler_init.h>
#include <pstl/execution>
#include <pstl/algorithm>
#else
#include <algorithm>
#endif
// External
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <map>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
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
  SORT(ret.begin(), ret.end());
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

  SORT(tfs_mapped.begin(), tfs_mapped.end());
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
          std::vector< std::pair< uint32_t, uint32_t> >& tneg,
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
    {
      tp++;
    }
    else
    {
      if (tneg.size() > 0)
      {
        if (std::binary_search(tneg.begin(), tneg.end(), x))
        {
          fp++;
        }
      }
      else
      {
        fp++;
      }
    }
  }
  return std::pair<uint32_t, uint32_t>(tp, fp);
}

void resize_by_fraction(std::vector< std::pair< uint32_t, uint32_t> >& truth,
                        std::vector< std::pair< uint32_t, uint32_t> >& tneg,
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
    {
      tp++;
    }
    else
    {
      if (tneg.size() > 0)
      {
        if (std::binary_search(tneg.begin(), tneg.end(), x))
        {
          fp++;
        }
      }
      else
      {
        fp++;
      }
    }

    cnt++;

    double f = numeric_cast<double> (tp) / numeric_cast<double> (truth.size());
    if (f >= fedg)
    {
      log(LOG_INFO) << "Resizing to top " << cnt << " egdes in order to keep "
                    << "at least " << truth.size() * f << " TP edges.\n";
      v.resize(cnt);
      break;
    }
  }

}

double print_roc(std::vector<std::pair<uint32_t, uint32_t>>& truth,
                 std::vector<std::pair<uint32_t, uint32_t>>& tneg,
                 std::vector<MiniEdge>& v,
                 std::pair<uint32_t, uint32_t>& tpfp,
                 uint32_t& datap,
                 std::ostream& out,
                 std::string algorithm)
{
  logger log(std::cerr, "ROC");
  uint32_t tp = 0, fp = 0;
  uint32_t cnt = 0;
  uint32_t ne = v.size();
  uint32_t intervx = 0;
  uint32_t intervy = 0;
  if (datap != 0)
  {
    intervy = numeric_cast<double> (ne) /
              numeric_cast<double> (datap);
    log(LOG_DEBUG) << "intervy: " << intervy << '\n';
  }
  out << "#TP/FP\t" <<  tpfp.first << '\t' << tpfp.second << '\n';

  // For on-the-fly integration -> AUC
  double y_i = 0, y_j = 0;
  double x_i = 0, x_j = 0;
  double auc = 0;

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
    {
      tp++;
    }
    else
    {
      if (tneg.size() > 0)
      {
        if (std::binary_search(tneg.begin(), tneg.end(), x))
        {
          fp++;
        }
      }
      else
      {
        fp++;
      }
    }
    double tpr = numeric_cast<double> (tp) / numeric_cast<double> (tpfp.first);
    double fpr = numeric_cast<double> (fp) / numeric_cast<double> (tpfp.second);
    double ppv = numeric_cast<double> (tp) / numeric_cast<double> (cnt);
    if (cnt > intervx || cnt == ne - 1)
    {
      out << tpr << '\t' << fpr << '\t' << ppv;
      if (algorithm != "")
      {
        out << '\t' << algorithm << '\n';
      }
      else
      {
        out << '\n';
      }
      intervx += intervy;
    }
    y_i = y_j;
    x_i = x_j;
    y_j = tpr;
    x_j = fpr;
    if (cnt > 1)
    {
      auc += ((y_i + y_j) / 2) * (x_j - x_i);
    }
  }
  return auc;
}

int roc(int argc, char * argv[])
{

  logger log(std::cerr, "ROC");

  // In case we have TBB/PSTL, initialize a global control object
  INIT_TBB_CONTROL();

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " roc";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  seidr_roc_param_t param;

  po::options_description
  umbrella("Calculate ROC curves of predictions in SeidrFiles given true edges");

  po::options_description opt("Common Options");
  opt.add_options()
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("help,h", "Show this help message")
  ("outfile,o",
   po::value<std::string>(&param.outfile)->default_value("-"),
   "Output file name ['-' for stdout]");

  po::options_description ompopt("OpenMP Options");
  ompopt.add_options()
  ("threads,O", po::value<int>(&param.nthreads)->
   default_value(GET_MAX_PSTL_THREADS()),
   "Number of OpenMP threads for parallel sorting");

  po::options_description ropt("ROC Options");
  ropt.add_options()
  ("index,i",
   po::value<uint32_t>(&param.tpos)->default_value(0, "last score"),
   "Index of score to use")
  ("edges,e",
   po::value<uint32_t>(&param.nedg)->default_value(0, "all"),
   "Number of top edges to consider")
  ("fraction,E",
   po::value<double>(&param.fedg)->default_value(-1, "all"),
   "Fraction of gold standard edges to include")
  ("points,p",
   po::value<uint32_t>(&param.datap)->default_value(0, "all"),
   "Number of data points to print")
  ("tfs,t",
   po::value<std::string>(&param.tfs),
   "List of transcription factors to consider")
  ("neg,x",
   po::value<std::string>(&param.gold_neg)->default_value("", ""),
   "True negative edges")
  ("all,a", po::bool_switch(&param.all)->default_value(false),
   "Calculate ROC for all scores in the SeidrFile")
  ;

  po::options_description req("Required");
  req.add_options()
  ("gold,g",
   po::value<std::string>(&param.gold),
   "Gold standard (true edges) input file")
  ("network,n", po::value<std::string>(&param.infile)->required(),
   "Input SeidrFile");

  umbrella.add(req).add(ropt).add(ompopt).add(opt);

  po::positional_options_description p;
  p.add("network", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, args).
            options(umbrella).positional(p).run(), vm);

  if (vm.count("help") || argc == 1)
  {
    std::cerr << umbrella << '\n';
    return 22;
  }

  try
  {
    po::notify(vm);
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  try
  {
    set_pstl_threads(param.nthreads, tbb_control);
    param.infile = to_absolute(param.infile);
    assert_exists(param.infile);
    assert_can_read(param.infile);

    param.gold = to_absolute(param.gold);
    assert_exists(param.gold);
    assert_can_read(param.gold);

    if (param.gold_neg != "")
    {
      param.gold_neg = to_absolute(param.gold_neg);
      assert_exists(param.gold_neg);
      assert_can_read(param.gold_neg);
    }

    if (param.outfile != "-")
    {
      param.outfile = to_absolute(param.outfile);
      if (! param.force)
      {
        assert_no_overwrite(param.outfile);
      }
      assert_dir_is_writeable(dirname(param.outfile));
    }
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  try
  {
    SeidrFile f(param.infile.c_str());
    f.open("r");
    SeidrFileHeader h;
    h.unserialize(f);
    make_tpos(param.tpos, h);
    std::map<std::string, uint32_t> refmap;
    uint32_t ind = 0;
    for (auto& no : h.nodes)
    {
      refmap[no] = ind++;
    }

    auto gv = make_gold_vec(param.gold, refmap);
    std::vector< std::pair< uint32_t, uint32_t > > ngv;
    if (param.gold_neg != "")
    {
      ngv = make_gold_vec(param.gold_neg, refmap);
    }

    std::shared_ptr<std::ostream> out;
    if (param.outfile == "-")
    {
      out = std::shared_ptr < std::ostream > (&std::cout, [](void*) {});
    }
    else
    {
      out = std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    // Store start of edges in the SF in case we need to read more than once
    int64_t offset = f.tell();
    if (! param.all)
    {
      // Start ROC for single algorithm
      auto nv = read_network_minimal(h, f, -1, param.tpos);
      SORTWCOMP(nv.begin(), nv.end(), me_score_sort_abs); //TODO: sort on rank, this is not robust
      if (param.tfs != "")
      {
        nv = filter_tf(param.tfs, nv, refmap);
      }
      if (param.nedg != 0)
      {
        nv.resize(param.nedg);
      }
      if (param.fedg > 0)
      {
        resize_by_fraction(gv, ngv, nv, param.fedg);
      }
      auto tpfp = get_tp_fp(gv, ngv, nv);
      double auc = print_roc(gv, ngv, nv, tpfp, param.datap, *out,
                             h.algs[param.tpos]);
      *out << "#AUC: " << auc << '\t' << h.algs[param.tpos] << '\n';
    }
    else
    {
      std::vector< std::pair<double, std::string> > aucs;
      for (uint16_t i = 0; i < h.attr.nalgs; i++)
      {
        f.seek(offset);
        auto nv = read_network_minimal(h, f, -1, i);
        SORTWCOMP(nv.begin(), nv.end(), me_score_sort_abs); //TODO: sort on rank, this is not robust
        if (param.tfs != "")
        {
          nv = filter_tf(param.tfs, nv, refmap);
        }
        if (param.nedg != 0)
        {
          nv.resize(param.nedg);
        }
        if (param.fedg > 0)
        {
          resize_by_fraction(gv, ngv, nv, param.fedg);
        }
        auto tpfp = get_tp_fp(gv, ngv, nv);
        double auc = print_roc(gv, ngv, nv, tpfp, param.datap, *out,
                               h.algs[i]);
        aucs.push_back( std::make_pair(auc, h.algs[i]) );
      }
      for (const auto& a : aucs)
      {
        *out << "#AUC: " << a.first << '\t' << a.second << '\n';
      }
    }

    f.close();
    return 0;
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

}
