#include <compare_clusters.h>
#include <common.h>
#include <fs.h>
#include <BSlogger.h>

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <tclap/CmdLine.h>

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

  std::string in_left;
  std::string in_right;
  std::string delim;

  double alpha = 0.05;
  unsigned long min_members = 20;
  unsigned long max_members = 100;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " cluster enrichment";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Test wether clusters of two networks show"
                       " significant overlap or extract clusters", ' ', 
                       version);

    TCLAP::ValueArg<std::string>
    arg_left("1", "first", "First cluster->gene mapping", true,
             "", "");
    cmd.add(arg_left);

    TCLAP::ValueArg<std::string>
    arg_right("2", "second", "Second cluster->gene mapping", false,
              "", "");
    cmd.add(arg_right);

    TCLAP::ValueArg<std::string>
    arg_delim("d", "delim", "Output delimiter", false,
              ",", ",");
    cmd.add(arg_delim);

    TCLAP::ValueArg<double>
    arg_alpha("a", "alpha", "Adjusted p-value cutoff", false,
              0.1, "0.1");
    cmd.add(arg_alpha);

    TCLAP::ValueArg<unsigned long>
    arg_minm("m", "min-members", "Minimum members of a cluster", false,
             20, "20");
    cmd.add(arg_minm);

    TCLAP::ValueArg<unsigned long>
    arg_maxm("M", "max-members", "Maximum members of a cluster", false,
             100, "100");
    cmd.add(arg_maxm);

    // Parse arguments
    cmd.parse(argc, args);
    in_left = arg_left.getValue();
    in_right = arg_right.getValue();
    alpha = arg_alpha.getValue();
    min_members = arg_minm.getValue();
    max_members = arg_maxm.getValue();
    delim = arg_delim.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return EINVAL;
  }

  if (in_right != "")
  {

    try
    {
      in_left = to_absolute(in_left);
      in_right = to_absolute(in_right);

      if (! file_can_read(in_left))
        throw std::runtime_error("Cannot read file " + in_left);

      if (! file_can_read(in_right))
        throw std::runtime_error("Cannot read file " + in_right);
    }
    catch (std::exception& e)
    {
      log(LOG_ERR) << e.what() << '\n';
      return 1;
    }

    auto x = cluster_map(in_left, min_members, max_members);
    auto y = cluster_map(in_right, min_members, max_members);

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

    std::cout << "First\tSecond\tFirst_Size\tSecond_size"
              << "\tMt\tNt\tMpat\tNpat\tRepr\tP_Adj\n";

    for (unsigned long i = 0; i < padj.size(); i++)
    {
      if (padj[i] < alpha)
      {
        std::cout << lname[i] << '\t'
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
    try
    {
      in_left = to_absolute(in_left);
      if (! file_can_read(in_left))
        throw std::runtime_error("Cannot read file " + in_left);
    }
    catch (std::exception& e)
    {
      log(LOG_ERR) << e.what() << '\n';
      return 1;
    }
    auto x = cluster_map(in_left, min_members, max_members);
    for (auto& i : x)
    {
      auto members = i.second.get_members();
      if (members.size() >= min_members && members.size() <= max_members)
      {
        std::cout << i.first << '\t';
        for (uint64_t j = 0; j < members.size(); j++)
        {
          std::cout << members[j]
                    << (j == members.size() - 1 ? "\n" : delim);
        }
      }
    }
  }
  return 0;
}

