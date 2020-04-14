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
#include <BSlogger.hpp>
#include <Serialize.h>
#include <common.h>
#include <dna.h>
#include <fs.h>
// External
#include <cerrno>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>

std::vector<std::string>
get_targets(const std::string& p)
{
  logger log(std::cerr, "dna");
  std::vector<std::string> ret;
  if (file_exists(p)) {
    // Get targets from file
    log(LOG_INFO) << "Getting targets from file\n";
    std::istream ifs(p.c_str(), std::ios::in);
    std::string line;
    while (std::getline(ifs, line)) {
      ret.push_back(line);
    }
  } else {
    // Get targets from command line
    log(LOG_INFO) << "Parsing targets from command line\n";
    ret = tokenize_delim(p, ',');
  }
  return ret;
}

std::unordered_map<std::string, uint64_t>
gene_index_map(const SeidrFileHeader& h)
{
  std::unordered_map<std::string, uint64_t> ret;
  uint64_t i = 0;
  for (const auto& g : h.nodes) {
    ret[g] = i++;
  }
  return ret;
}

std::vector<std::string>
gene_intersect(const std::vector<std::string>& v1,
               const std::vector<std::string>& v2,
               bool sorted = false)
{
  std::vector<std::string> ret;
  if (!sorted) {
    std::vector v1_copy = v1;
    SORT(v1_copy.begin(), v1_copy.end());
    for (const auto& g : v2) {
      if (std::binary_search(v1_copy.begin(), v1_copy.end(), g)) {
        ret.push_back(g);
      }
    }
  } else {
    for (const auto& g : v2) {
      if (std::binary_search(v1.begin(), v1.end(), g)) {
        ret.push_back(g);
      }
    }
  }
  return ret;
}

std::unordered_map<std::string, std::vector<string>>
pathway_members(const std::string& gene_to_pw)
{
  std::unordered_map<std::string, std::vector<string>> ret;
  std::ifstream ifs(gene_to_pw.c_str());
  std::string line;
  while (std::getline(ifs, line)) {
    std::string gene;
    std::string pathways;
    std::stringstream ss(line);
    ss >> gene;
    ss >> pathways;
    std::vector<std::string> pids = tokenize_delim(pathways, '|');
    for (const auto& pid : pids) {
      ret[pid].push_back(gene);
    }
  }
  return ret;
}

double
diff_score(
  const seidr_dna_param_t& param,
  const std::string& pathway,
  const std::unordered_map<std::string, std::vector<string>>& pathway_members,
  const std::unordered_map<std::string, uint64_t>& gene_index1,
  const std::unordered_map<std::string, uint64_t>& gene_index2,
  const arma::mat& X1,
  const arma::mat& X2,
  const std::vector<std::string>& sorted_isect)
{
  // TODO try catch this
  const std::vector<std::string>& members = pathway_members.at(pathway);
  std::vector<std::string> filtered_members;
  for (const auto& g : members) {
    if (std::binary_search(sorted_isect.begin(), sorted_isect.end(), g)) {
      filtered_members.push_back(g);
    }
  }
  // TODO get d0

  // Get N permutations
}

int
dna(int argc, char* argv[])
{

  logger log(std::cerr, "dna");

  // Variables used by the function
  seidr_dna_param_t param;

  // We ignore the first argument
  const char* args[argc - 1];
  std::string pr(argv[0]);
  pr += " dna";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++)
    args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella("Perform pathway mediated differential"
                                   " network analysis.");

  po::options_description opt("Common Options");
  opt.add_options()("force,f",
                    po::bool_switch(&param.force)->default_value(false),
                    "Force overwrite output file if it exists")(
    "help,h", "Show this help message")(
    "outfile,o",
    po::value<std::string>(&param.file_out)->default_value("-"),
    "Output file name ['-' for stdout]");

  po::options_description copt("DNA Options");
  copt.add_options()(
    "index-a,i",
    po::value<uint32_t>(&param.tindex1)->default_value(0, "last score"),
    "Use scores on this index for network A")(
    "index-b,j",
    po::value<uint32_t>(&param.tindex1)->default_value(0, "last score"),
    "Use scores on this index for network B")(
    "n-permutations,B",
    po::value<uint64_t>(&param.nperm)->default_value(1000),
    "The number of permutations to use when calculating the differential "
    "score"),
    ("target,t",
     po::value<std::string>(&param.target),
     "Calculate P-value for these pathways (can be multiple comma separated or "
     "a "
     "file with one target per line)");

  po::options_description req("Required Options [can be positional]");
  req.add_options()("gene-to-pathway",
                    po::value<std::string>(&param.gene_to_pw)->required(),
                    "Gene to pathway mapping.")(
    "network-1",
    po::value<std::string>(&param.net1)->required(),
    "Input SeidrFile for network A")(
    "network-2",
    po::value<std::string>(&param.net2)->required(),
    "Input SeidrFile for network B");

  umbrella.add(req).add(copt).add(opt);

  po::positional_options_description p;
  p.add("gene-to-pathway", 1);
  p.add("network-1", 1);
  p.add("network-2", 1);

  po::variables_map vm;
  po::store(
    po::command_line_parser(argc, args).options(umbrella).positional(p).run(),
    vm);

  if (vm.count("help") || argc == 1) {
    std::cerr << umbrella << '\n';
    return EINVAL;
  }

  try {
    po::notify(vm);
  } catch (std::exception& e) {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }

  param.net1 = to_absolute(param.net1);
  param.net2 = to_absolute(param.net2);

  try {
    assert_exists(param.net1);
    assert_exists(param.net2);
    assert_exists(param.gene_to_pw);
    assert_can_read(param.net1);
    assert_can_read(param.net2);
    assert_can_read(param.gene_to_pw);
    if (!param.force && param.file_out != "-") {
      assert_no_overwrite(param.file_out);
      assert_dir_is_writeable(dirname(param.file_out));
    }
  } catch (std::runtime_error& e) {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

  SeidrFile n1(param.net1.c_str());
  n1.open("r");
  SeidrFile n2(param.net2.c_str());
  n2.open("r");
  SeidrFileHeader h1;
  h1.unserialize(n1);
  SeidrFileHeader h2;
  h2.unserialize(n2);
  make_tpos(param.tindex1, h1);
  make_tpos(param.tindex2, h2);

  log(LOG_INFO) << "Using scores from " << h1.algs[param.tindex1]
                << " for network 1 and " << h2.algs[param.tindex2]
                << " for network 2.\n";

  log(LOG_INFO) << "Reading network 1 into memory\n";
  arma::mat X1 = read_network_arma(
    h1, n1, -std::numeric_limits<double>::infinity(), param.tindex1, false);
  log(LOG_INFO) << "Reading network 2 into memory\n";
  arma::mat X2 = read_network_arma(
    h2, n2, -std::numeric_limits<double>::infinity(), param.tindex2, false);

  std::vector<std::string> sorted_isect =
    gene_intersect(h1.nodes(), h2.nodes());
  SORT(sorted_isect.begin(), sorted_isect.end());

  n1.close();
  n2.close();
  return 0;
}