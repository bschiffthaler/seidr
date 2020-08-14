/* * Seidr - Create and operate on gene crowd networks
 * Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
 *
 * This file is part of Seidr.
 *
 * Seidr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Seidr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

constexpr double SEIDR_COMPARE_CLUST_DEF_ALPHA = 0.1;
constexpr uint64_t SEIDR_COMPARE_CLUST_DEF_MIN_MEMBERS = 20;
constexpr uint64_t SEIDR_COMPARE_CLUST_DEF_MAX_MEMBERS = 200;

struct seidr_cc_param_t
{
  std::string in_left;
  std::string in_right;
  std::string delim;
  double alpha = SEIDR_COMPARE_CLUST_DEF_ALPHA;
  uint64_t min_members = SEIDR_COMPARE_CLUST_DEF_MIN_MEMBERS;
  uint64_t max_members = SEIDR_COMPARE_CLUST_DEF_MAX_MEMBERS;
  std::string file_out;
  bool force;
};

struct result
{
  double pval = 0;
  uint64_t mpat = 0;
  uint64_t npat = 0;
  uint64_t mt = 0;
  uint64_t nt = 0;
  uint64_t lhs = 0;
  uint64_t rhs = 0;
};

class cluster
{
public:
  // getters
  std::string get_name() { return _name; }
  std::vector<std::string>& get_members() { return _members; }
  std::string get_parent_name() { return _parent_name; }
  std::vector<std::string> parent_members();
  std::shared_ptr<cluster> get_parent() { return _parent; }
  // setters
  void set_name(std::string& n) { _name = n; }
  void add_member(std::string& m) { _members.push_back(m); }
  void set_parent_name(std::string& c)
  {
    _parent_name = c;
    _has_parent = true;
  }
  void set_parent(cluster& c)
  {
    _parent = std::make_shared<cluster>(c);
    _has_parent = true;
  }
  // logical
  bool has_parent() const { return _has_parent; }
  // utility
  void sort_members();
  void make_unique();

private:
  std::string _name;
  bool _has_parent;
  std::vector<std::string> _members;
  std::string _parent_name;
  std::shared_ptr<cluster> _parent;
};

inline void
cluster::sort_members()
{
  std::sort(_members.begin(), _members.end());
}

inline void
cluster::make_unique()
{
  auto it = std::unique(_members.begin(), _members.end());
  _members.resize(std::distance(_members.begin(), it));
}

struct prank
{
  uint64_t orig_pos;
  double cmin;
  double val;
  friend bool operator<(const prank& a, const prank& b);
};

inline bool
operator<(const prank& a, const prank& b)
{
  return a.val > b.val;
}

double
log10_sigma_t(uint64_t m, uint64_t n, uint64_t mt, uint64_t nt);
result
test(cluster& lhs, cluster& rhs);
double
log10_binom_coef(uint64_t n, uint64_t k);
std::vector<double>
bh_adjust(std::vector<double> pvals);
int
cluster_enrichment(const std::vector<std::string>& args);
double
div_as_double(uint64_t lhs, uint64_t rhs);
