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

#include <compare_clusters.h>
#include <common.h>
#include <fs.h>
#include <BSlogger.hpp>

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

double div_as_double(unsigned long lhs,
                     unsigned long rhs)
{
  double a = lhs;
  double b = rhs;
  return a / b;
}

std::vector<std::string> strsplit(std::string& sin, char split)
{
  std::vector<std::string> ret;
  std::string word{};
  for (auto& c : sin)
  {
    if (c == split)
    {
      if (word.size() > 0)
        ret.push_back(word);
    }
    word += c;
  }
  ret.push_back(word);
  return ret;
}

std::map<std::string, cluster>
cluster_map(std::string inf,
            unsigned long min_members,
            unsigned long max_members)
{
  std::map<std::string, cluster> ret;
  std::map<std::string, cluster> filt;
  std::ifstream ifs(inf.c_str());
  std::string line;

  // Read file and add clusters to map
  while (std::getline(ifs, line))
  {
    std::string field;
    std::stringstream ss(line);
    ss >> field;
    std::vector<std::string> clusters = strsplit(field, ':');
    ss >> field;
    for (unsigned long i = 0; i < clusters.size(); i++)
    {
      auto hit = ret.find(clusters[i]);
      if (hit == ret.end())
      {
        ret[clusters[i]] = cluster();
        ret[clusters[i]].set_name(clusters[i]);
        if (i > 0)
        {
          ret[clusters[i]].set_parent_name(clusters[i - 1]);
        }
      }
      ret[clusters[i]].add_member(field);
    }
  }

  // Drop modules with too few/many members
  for (auto it = ret.begin(); it != ret.end(); it++)
  {
    if ((it->second.get_members().size() >= min_members) &&
        (it->second.get_members().size() <= max_members))
    {
      //std::cout << it->first << ' ' << it->second.get_members().size() << '\n';
      filt[it->first] = it->second;
    }
  }

  // Ensure vectors are sorted and unique
  for (auto it = filt.begin(); it != filt.end(); it++)
  {
    it->second.sort_members();
    it->second.make_unique();
  }

  // Link pointers to parents
  for (auto it = filt.begin(); it != filt.end(); it++)
  {
    if (it->second.has_parent())
    {
      it->second.set_parent( filt[it->second.get_parent_name()] );
    }
  }
  return filt;
}

std::vector<std::string> cluster::parent_members()
{
  std::vector<std::string> members = this->get_members();
  auto current_cluster = (*this);
  while (current_cluster.has_parent())
  {
    auto parent = current_cluster.get_parent();
    current_cluster = *parent;
    for (auto& m : current_cluster.get_members())
      members.push_back(m);
  }
  std::sort(members.begin(), members.end());
  auto it = std::unique(members.begin(), members.end());
  members.resize(std::distance(members.begin(), it));
  return members;
}

double log10_binom_coef(unsigned long n, unsigned long k)
{
  if (k > n) {
    return -1;
  }
  if (k == 0) return 0;

  double nd = n;

  double r = log10( nd );

  for (unsigned long i = 2; i <= k; i++)
  {
    double id =  i;
    r += log10( ( nd + 1 - id ) / id );
  }
  return r;
}

unsigned long min(unsigned long a, unsigned long b)
{
  return a <= b ? a : b;
}

// Probability to observe sigma_t or more terms, equ (2) in
// Grossman et al.
double log10_sigma_t(unsigned long m, unsigned long n,
                     unsigned long mt, unsigned long nt)
{
  unsigned long limit = min(mt, n);

  double s = 0;

  for (unsigned long k = nt; k <= limit; k++)
  {
    // We get logarithms from log10_binom_coef so products
    // become sums and divisions differences
    double t1 = log10_binom_coef(mt, k);
    double t2 = log10_binom_coef(m - mt, n - k);
    double t3 = log10_binom_coef(m, n);

    if (t1 < 0 || t2 < 0) continue;
    if (t3 < 0) return -1;

    double e = t1 + t2 - t3;
    s += pow(10, e);
  }

  return s;
}

result test(cluster& lhs, cluster& rhs)
{
  unsigned long mt = rhs.get_members().size();

  std::vector<std::string> lhs_members = lhs.get_members();
  std::vector<std::string> rhs_members = rhs.get_members();
  std::vector<std::string> rhs_parents = rhs.parent_members();

  unsigned long npat = 0;
  unsigned long nt = 0;

  for (auto& n : lhs_members)
    if (std::binary_search(rhs_parents.begin(), rhs_parents.end(), n))
      npat++;
  for (auto& n : lhs_members)
    if (std::binary_search(rhs_members.begin(), rhs_members.end(), n))
      nt++;

  unsigned long mpat = rhs_parents.size();

  //std::cerr << lhs.get_members().size() << '(' << npat << "), "
  //          << mt << '(' << mpat << "), " << nt << "\t";

  result res;
  res.npat = npat;
  res.mpat = mpat;
  res.mt = mt;
  res.nt = nt;
  res.pval = log10_sigma_t(mpat, npat, mt, nt);
  res.lhs = lhs_members.size();
  res.rhs = rhs_members.size();
  return res;
}

std::vector<double> bh_adjust(std::vector< double > pvals)
{
  std::vector< prank > pr(pvals.size());
  std::vector< double > ret(pvals.size());
  //Assign p values and sort descending
  for ( unsigned long i = 0; i < pvals.size(); i++)
  {
    pr[i].orig_pos = i;
    pr[i].val = pvals[i];
  }
  std::sort(pr.begin(), pr.end());

  //Process first entry
  double n = pvals.size();
  pr[0].cmin = pr[0].val * 1;
  if (pr[0].cmin > 1) pr[0].cmin = 1;
  ret[pr[0].orig_pos] = pr[0].cmin;

  double prev = pr[0].cmin;

  //Process rest
  for (unsigned long i = 1; i < pvals.size(); i++)
  {
    double ix = i;
    double cur = n / (n - ix) * pr[i].val;
    if (cur < prev)
    {
      pr[i].cmin = cur > 1 ? 1 : cur;
      prev = cur;
    }
    else
    {
      pr[i].cmin = prev;
    }

    ret[pr[i].orig_pos] = pr[i].cmin;
  }
  return ret;
}

int cluster_enrichment(int argc, char ** argv)
{

  logger log(std::cerr, "cluster_enrichment");

  seidr_cc_param_t param;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " cluster enrichment";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella("Test wether clusters of two networks show"
                                   " significant overlap or extract clusters");

  po::options_description opt("Common Options");
  opt.add_options()
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("help,h", "Show this help message")
  ("outfile,o", po::value<std::string>(&param.file_out)->default_value("-"),
   "Output file name ['-' for stdout]");

  po::options_description ccopt("Cluster-Compare Options");
  ccopt.add_options()
  ("second,2",
   po::value<std::string>(&param.in_right),
   "Another cluster->gene mapping [can be positional]")
  ("delim,d",
   po::value<std::string>(&param.delim)->default_value(","),
   "Output delimiter")
  ("alpha,a",
   po::value<double>(&param.alpha)->default_value(0.1, "0.1"),
   "Adjusted p-value cutoff")
  ("min-members,m",
   po::value<unsigned long>(&param.min_members)->default_value(20),
   "Minimum members of a cluster")
  ("max-members,M",
   po::value<unsigned long>(&param.max_members)->default_value(200),
   "Maximum members of a cluster");

  po::options_description req("Required Options [can be positional]");
  req.add_options()
  ("first,1",
   po::value<std::string>(&param.in_left)->required(),
   "First cluster->gene mapping");

  umbrella.add(req).add(ccopt).add(opt);

  po::positional_options_description p;
  p.add("first", 1);
  p.add("second", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, args).
            options(umbrella).positional(p).run(), vm);

  if (vm.count("help") || argc == 1)
  {
    std::cerr << umbrella << '\n';
    return EINVAL;
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
    param.in_left = to_absolute(param.in_left);
    assert_exists(param.in_left);
    assert_can_read(param.in_left);
    if (param.file_out != "-")
    {
      param.file_out = to_absolute(param.file_out);
      assert_dir_is_writeable(dirname(param.file_out));
      if (! param.force)
      {
        assert_no_overwrite(param.file_out);
      }
    }
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  std::shared_ptr<std::ostream> out;

  if (param.file_out == "-")
  {
    out.reset(&std::cout, [](...) {});
  }
  else
  {
    out.reset(new std::ofstream(param.file_out));
  }

  if (vm.count("second"))
  {
    try
    {
      param.in_right = to_absolute(param.in_right);
      assert_exists(param.in_right);
      assert_can_read(param.in_right);
    }
    catch (std::exception& e)
    {
      log(LOG_ERR) << e.what() << '\n';
      return 1;
    }

    auto x = cluster_map(param.in_left, param.min_members, param.max_members);
    auto y = cluster_map(param.in_right, param.min_members, param.max_members);

    std::vector<std::string> lname, rname;
    std::vector<double> padj;
    std::vector<unsigned long> mpat, npat, mt, nt, lhss, rhss;

    for (auto& lhs : x)
    {
      for (auto& rhs : y)
      {
        lname.push_back(lhs.first);
        rname.push_back(rhs.first);
        result res = test(lhs.second, rhs.second);
        padj.push_back(res.pval);
        mpat.push_back(res.mpat);
        npat.push_back(res.npat);
        mt.push_back(res.mt);
        nt.push_back(res.nt);
        lhss.push_back(res.lhs);
        rhss.push_back(res.rhs);
      }
    }

    padj = bh_adjust(padj);

    (*out) << "First\tSecond\tFirst_Size\tSecond_size"
           << "\tMt\tNt\tMpat\tNpat\tRepr\tP_Adj\n";

    for (unsigned long i = 0; i < padj.size(); i++)
    {
      if (padj[i] < param.alpha)
      {
        (*out) << lname[i] << '\t'
               << rname[i] << '\t'
               << lhss[i] << '\t'
               << rhss[i] << '\t'
               << mt[i] << '\t'
               << nt[i] << '\t'
               << mpat[i] << '\t'
               << npat[i] << '\t'
               << div_as_double(nt[i], rhss[i]) << '\t'
               << padj[i] << '\n';
      }
    }
  }
  else
  {
    auto x = cluster_map(param.in_left, param.min_members, param.max_members);
    for (auto& i : x)
    {
      auto members = i.second.get_members();
      if (members.size() >= param.min_members &&
          members.size() <= param.max_members)
      {
        (*out) << i.first << '\t';
        for (uint64_t j = 0; j < members.size(); j++)
        {
          (*out) << members[j]
                 << (j == members.size() - 1 ? "\n" : param.delim);
        }
      }
    }
  }
  return 0;
}

