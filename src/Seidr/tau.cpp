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

// Calculate pairwise tau scores for nodes in two networks
// Mostly from SciPy https://github.com/scipy/scipy/blob/v0.15.1/scipy/stats/stats.py#L2827
// and https://www.geeksforgeeks.org/merge-sort/

// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
#include <tau.h>
// External
#include <armadillo>
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

uint64_t ties(const arma::vec& arr)
{
  uint64_t ret = 0;
  uint64_t l = 0;

  for (uint64_t i = 1; i < arr.n_elem; i++)
  {
    if (! (almost_equal(arr(l), arr(i))))
    {
      ret += ((i - l) * (i - l - 1)) / 2;
      l = i;
    }
  }

  ret += ((arr.n_elem - l) * (arr.n_elem - l - 1)) / 2;

  return ret;
}

uint64_t joint_ties(const std::vector<std::pair<double, double>>& arr)
{

  uint64_t ret = 0;
  uint64_t l = 0;

  for (uint64_t i = 1; i < arr.size(); i++)
  {
    if ((( ! almost_equal(arr[l].first, arr[i].first) ) ||
         ( ! almost_equal(arr[l].second, arr[i].second))))
  {
    ret += ((i - l) * (i - l - 1)) / 2;
      l = i;
    }
  }

  ret += ((arr.size() - l) * (arr.size() - l - 1)) / 2;

  return ret;
}

bool comp_pd(const std::pair<double, double>& lhs,
             const std::pair<double, double>& rhs)
{
  if (almost_equal(lhs.first, rhs.first))
  {
    return lhs.second < rhs.second;
  }
  return lhs.first < rhs.first;
}

uint64_t merge(arma::vec& arr, arma::uword l, arma::uword m, arma::uword r)
{

  uint64_t swaps = 0;

  arma::uword n1 = m - l + 1;
  arma::uword n2 = r - m;

  arma::vec L(n1);
  arma::vec R(n2);

  for (arma::uword i = 0; i < n1; i++)
    L(i) = arr(l + i);
  for (arma::uword j = 0; j < n2; j++)
    R(j) = arr(m + 1 + j);

  arma::uword i = 0;
  arma::uword j = 0;
  arma::uword k = l;

  while (i < n1 && j < n2)
  {
    if (L(i) <= R(j))
    {
      arr(k) = L(i);
      i++;
    }
    else
    {
      arr(k) = R(j);
      j++;
      swaps = swaps + m - i + 1 - l;
    }
    k++;
  }
  while (i < n1)
  {
    arr(k) = L(i);
    i++;
    k++;
  }
  while (j < n2)
  {
    arr(k) = R(j);
    j++;
    k++;
  }
  return swaps;
}

uint64_t  mergeSort(arma::vec& arr, arma::uword l, arma::uword r)
{
  uint64_t swaps = 0;
  if (l < r)
  {
    int m = l + (r - l) / 2;
    swaps = mergeSort(arr, l, m);
    swaps += mergeSort(arr, m + 1, r);
    swaps += merge(arr, l, m, r);
  }
  return swaps;
}

double k_tau_ms(const arma::vec& x, const arma::vec& y)
{
  std::vector<std::pair<double, double>> pairs;
  for (arma::uword i = 0; i < x.n_elem; i++)
  {
    pairs.push_back(std::pair<double, double>(x(i), y(i)));
  }
  // First sort by x, sort ties in x by y, using almost_equals for first pair
  std::sort(pairs.begin(), pairs.end(), comp_pd);
  // Now mergesort y only, keeping track of the number of swaps
  arma::vec xtmp(pairs.size());
  arma::vec ytmp(pairs.size());
  for (uint64_t i = 0; i < pairs.size(); i++)
  {
    xtmp(i) = pairs[i].first;
    ytmp(i) = pairs[i].second;
  }

  uint64_t t = joint_ties(pairs);
  uint64_t u = ties(xtmp);

  uint64_t swaps = mergeSort(ytmp, 0, ytmp.n_elem - 1);

  uint64_t v = ties(ytmp);
  uint64_t ltri = (x.n_elem * (x.n_elem - 1)) / 2;

  if (u == ltri || v == ltri)
    return std::numeric_limits<double>::quiet_NaN();

  double ltrid = ltri;
  double ud = u;
  double vd = v;
  double td = t;
  double swapsd = swaps;


  double denom = exp(0.5 * (log(ltrid - ud) + log(ltrid - vd)));
  double tau = ((ltrid - (vd + ud - td)) - 2.0 * swapsd) / denom;

  return tau;
}

void get_tau(const seidr_param_tau_t& param)
{
  uint32_t si1 = param.tpos_a;
  uint32_t si2 = param.tpos_b;

  std::shared_ptr<std::ostream> out;
  if (param.out_file == "-")
    out = std::shared_ptr < std::ostream > (&std::cout, [](void*) {});
  else
    out = std::shared_ptr<std::ostream>
          (new std::ofstream(param.out_file.c_str()));

  SeidrFile sf1(param.network_a.c_str());
  SeidrFile sf2(param.network_b.c_str());
  sf1.open("r");
  sf2.open("r");
  SeidrFileHeader h1;
  SeidrFileHeader h2;
  h1.unserialize(sf1);
  h2.unserialize(sf2);

  make_tpos(si1, h1);
  make_tpos(si2, h2);

  std::map<std::string, uint32_t> node_union;

  for (auto node : h1.nodes)
  {
    node_union[node] = 0;
  }
  for (auto node : h2.nodes)
  {
    node_union[node] = 0;
  }

  std::vector<std::string> colheader;
  for (auto n : node_union)
    colheader.push_back(n.first);

  uint32_t ctr = 0;
  for (auto& node : node_union)
  {
    node.second = ctr++;
  }

  arma::mat rm1(node_union.size(), node_union.size(), arma::fill::zeros);
  arma::mat rm2(node_union.size(), node_union.size(), arma::fill::zeros);

  sf1.seek(0);
  sf2.seek(0);

  sf1.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    std::string& i = h.nodes[e.index.i];
    std::string& j = h.nodes[e.index.j];
    double v = e.scores[si1].s;
    rm1(node_union[i], node_union[j]) = v;
    rm1(node_union[j], node_union[i]) = v;
  });

  sf2.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    std::string& i = h.nodes[e.index.i];
    std::string& j = h.nodes[e.index.j];
    double v = e.scores[si2].s;
    rm2(node_union[i], node_union[j]) = v;
    rm2(node_union[j], node_union[i]) = v;
  });

  sf1.close();
  sf2.close();

  for (arma::uword i = 0; i < node_union.size(); i++)
  {
  	arma::vec x(rm1.n_rows - 1);
  	arma::vec y(rm2.n_rows - 1);
  	arma::uword idx = 0;
  	for (arma::uword j = 0; j < rm1.n_rows; j++)
  	{
  		if (i != j)
  		{
  			x(idx) = rm1(j,i);
  			y(idx) = rm2(j,i);
  			idx++;
  		}
  	}
  	double tau = k_tau_ms(x, y);
  	std::cout << colheader[i] << '\t' << tau << '\n';
  }
}

void get_tau_cs(const seidr_param_tau_t& param)
{
  uint32_t si1 = param.tpos_a;
  uint32_t si2 = param.tpos_b;

  std::shared_ptr<std::ostream> out;
  if (param.out_file == "-")
    out = std::shared_ptr < std::ostream > (&std::cout, [](void*) {});
  else
    out = std::shared_ptr<std::ostream>
          (new std::ofstream(param.out_file.c_str()));

  SeidrFile sf1(param.network_a.c_str());
  SeidrFile sf2(param.network_b.c_str());
  sf1.open("r");
  sf2.open("r");
  SeidrFileHeader h1;
  SeidrFileHeader h2;
  h1.unserialize(sf1);
  h2.unserialize(sf2);

  make_tpos(si1, h1);
  make_tpos(si2, h2);

  std::map<std::string, std::set<std::string>> a_to_b;
  std::map<std::string, std::set<std::string>> b_to_a;

  std::ifstream dict_file(param.dict.c_str());
  std::string line;
  while (std::getline(dict_file, line))
  {
  	std::stringstream ss(line);
  	std::string lhs;
  	std::string rhs;
  	ss >> lhs;
  	ss >> rhs;
  	a_to_b[lhs].insert(rhs);
  	b_to_a[rhs].insert(lhs);
  }

  std::map<std::string, uint32_t> node_union;

  for (auto node : h1.nodes)
  {
  	// No orthologs
  	if (a_to_b.find(node) == a_to_b.end())
    	node_union[node] = 0;
  }
  for (auto node : h2.nodes)
  {
  	// No orthologs
  	if (b_to_a.find(node) == b_to_a.end())
    	node_union[node] = 0;
  }
  // Add orthologs 
  for (auto& ortho_source : a_to_b)
  {
  	for (auto& ortho_target : ortho_source.second)
  	{
  		node_union[ ortho_source.first + "," + ortho_target ] = 0;
  	}
  }
  // Add orthologs from second net inversing the order of IDs
  for (auto& ortho_source : b_to_a)
  {
  	for (auto& ortho_target : ortho_source.second)
  	{
  		node_union[ ortho_target + "," + ortho_source.first ] = 0;
  	}
  }

  std::vector<std::string> colheader;
  for (auto n : node_union)
    colheader.push_back(n.first);

  uint32_t ctr = 0;
  for (auto& node : node_union)
    node.second = ctr++;

  arma::mat rm1(node_union.size(), node_union.size(), arma::fill::zeros);
  arma::mat rm2(node_union.size(), node_union.size(), arma::fill::zeros);

  sf1.seek(0);
  sf2.seek(0);

  sf1.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    std::string i = h.nodes[e.index.i];
    std::string j = h.nodes[e.index.j];
    double v = e.scores[si1].s;

    auto ptr_i = a_to_b.find(i);
		auto ptr_j = a_to_b.find(j);

    if (ptr_i == a_to_b.end() && ptr_j == a_to_b.end())
    {
    	rm1(node_union[i], node_union[j]) = v;
    	rm1(node_union[j], node_union[i]) = v;
  	}
  	else if (ptr_i != a_to_b.end() && ptr_j != a_to_b.end()) 
  	{
  		for (auto si : ptr_i->second)
  		{
  			for (auto sj : ptr_j -> second)
  			{
  				auto xi = ptr_i->first + "," + si;
  				auto xj = ptr_j->first + "," + sj;
  				rm1(node_union.at(xi), node_union.at(xj)) = v;
    			rm1(node_union.at(xj), node_union.at(xi)) = v;
  			}
  		}
  	}
  	else if (ptr_i != a_to_b.end())
  	{
  		for (auto si : ptr_i->second)
  		{
  			auto xi = ptr_i->first + "," + si;
  			rm1(node_union.at(xi), node_union.at(j)) = v;
    		rm1(node_union.at(j), node_union.at(xi)) = v;
			}
  	}
  	else if (ptr_j != a_to_b.end())
  	{
  		for (auto sj : ptr_j->second)
  		{
  			auto xj = ptr_j->first + "," + sj;
  			rm1(node_union.at(i), node_union.at(xj)) = v;
    		rm1(node_union.at(xj), node_union.at(i)) = v;
			}
  	}
  });

  sf2.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    std::string& i = h.nodes[e.index.i];
    std::string& j = h.nodes[e.index.j];
    double v = e.scores[si2].s;

    auto ptr_i = b_to_a.find(i);
		auto ptr_j = b_to_a.find(j);

    if (ptr_i == a_to_b.end() && ptr_j == a_to_b.end())
    {
    	rm2(node_union[i], node_union[j]) = v;
    	rm2(node_union[j], node_union[i]) = v;
  	}
  	else if (ptr_i != a_to_b.end() && ptr_j != a_to_b.end()) 
  	{
  		for (auto si : ptr_i->second)
  		{
  			for (auto sj : ptr_j -> second)
  			{
  				auto xi = si + "," + ptr_i->first;
  				auto xj = sj + "," + ptr_j->first;
  				rm2(node_union.at(xi), node_union.at(xj)) = v;
    			rm2(node_union.at(xj), node_union.at(xi)) = v;
  			}
  		}
  	}
  	else if (ptr_i != a_to_b.end())
  	{
  		for (auto si : ptr_i->second)
  		{
  			auto xi = si + "," + ptr_i->first;
  			rm2(node_union.at(xi), node_union.at(j)) = v;
    		rm2(node_union.at(j), node_union.at(xi)) = v;
			}
  	}
  	else if (ptr_j != a_to_b.end())
  	{
  		for (auto sj : ptr_j->second)
  		{
  			auto xj = sj + "," + ptr_j->first;
  			rm2(node_union.at(i), node_union.at(xj)) = v;
    		rm2(node_union.at(xj), node_union.at(i)) = v;
			}
  	}
  });

  sf1.close();
  sf2.close();

  for (arma::uword i = 0; i < node_union.size(); i++)
  {
  	arma::vec x(rm1.n_rows - 1);
  	arma::vec y(rm2.n_rows - 1);
  	arma::uword idx = 0;
  	for (arma::uword j = 0; j < rm1.n_rows; j++)
  	{
  		if (i != j)
  		{
  			x(idx) = rm1(j,i);
  			y(idx) = rm2(j,i);
  			idx++;
  		}
  	}
  	double tau = k_tau_ms(x, y);
  	std::cout << colheader[i] << '\t' << tau << '\n';
  }
}

int tau(int argc, char * argv[]) {

  LOG_INIT_CERR();

  // Variables used by the function
  seidr_param_tau_t param;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " tau";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella;

  po::options_description opt("Common Options");
  opt.add_options()
  ("help,h", "Show this help message")
  ("force,f", po::bool_switch(&param.force)->default_value(false),
   "Force overwrite output file if it exists")
  ("out-file,o", po::value<std::string>(&param.out_file)->default_value("-"),
   "Output file name ['-' for stdout]");

  po::options_description topt("Tau Options");
  topt.add_options()
  ("dict,d", po::value<std::string>(&param.dict)->default_value(""),
   "Dictionary of orthologies")
  ("index-a,a", po::value<uint32_t>(&param.tpos_a)->default_value(0),
   "Score index to use for first network")
  ("index-b,b", po::value<uint32_t>(&param.tpos_b)->default_value(0),
   "Score index to use for second network");

  po::options_description req("Required");
  req.add_options()
  ("net-a", po::value<std::string>(&param.network_a)->required(),
   "First input SeidrFile")
  ("net-b", po::value<std::string>(&param.network_b)->required(),
   "Second input SeidrFile");

  umbrella.add(req).add(topt).add(opt);

  po::positional_options_description p;
  p.add("net-a", 1);
  p.add("net-b", 1);

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
    param.network_a = to_absolute(param.network_a);
    assert_exists(param.network_a);
    assert_can_read(param.network_a);

    param.network_b = to_absolute(param.network_b);
    assert_exists(param.network_b);
    assert_can_read(param.network_b);

    if (param.out_file != "-")
    {
      param.out_file = to_absolute(param.out_file);
      if (! param.force)
      {
        assert_no_overwrite(param.out_file);
      }
      assert_dir_is_writeable(dirname(param.out_file));
    }
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  if (param.dict == "")
  	get_tau(param);
  else
  	get_tau_cs(param);

  return 0;
}