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

struct seidr_cc_param_t
{
  std::string in_left;
  std::string in_right;
  std::string delim;
  double alpha = 0.05;
  unsigned long min_members = 20;
  unsigned long max_members = 100;
  std::string file_out;
  bool force;
};

struct result
{
  double pval = 0;
  unsigned long mpat = 0;
  unsigned long npat = 0;
  unsigned long mt = 0;
  unsigned long nt = 0;
  unsigned long lhs = 0;
  unsigned long rhs = 0;
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
  bool has_parent() { return _has_parent; }
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
  unsigned long orig_pos;
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
log10_sigma_t(unsigned long m,
              unsigned long n,
              unsigned long mt,
              unsigned long nt);
result
test(cluster& lhs, cluster& rhs);
double
log10_binom_coef(unsigned long n, unsigned long k);
std::vector<double>
bh_adjust(std::vector<double> pvals);
int
cluster_enrichment(int argc, char** argv);
double
div_as_double(unsigned long lhs, unsigned long rhs);
