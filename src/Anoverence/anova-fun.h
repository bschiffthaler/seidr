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
#include <armadillo>
#include <set>
#include <string>
#include <vector>

typedef arma::uword seidr_index;

#define ANOVA_FULL 0
#define ANOVA_PARTIAL 1

struct seidr_anova_param_t
{
  std::string infile;
  float weight;
  std::string feature_file;
  size_t mode = ANOVA_FULL;
  std::string outfile;
  bool force = false;
  std::string targets_file;
  std::string gene_file;
  unsigned verbosity;
};

struct sparse
{

  std::string index; // perturbation ID
  float value = 0;   // weight of the perturbation

  sparse(std::string idx, float val)
    : index(idx)
    , value(val)
  {}

  bool operator<(const sparse& s2) const
  {
    if (index < s2.index) {
      return true;
    } else if (index == s2.index) {
      if (value < s2.value)
        return true;
    }
    return false;
  }
  bool operator==(const sparse& s2) const
  {
    return (index == s2.index) && (value == s2.value);
  }
};

struct gene
{
  int over = 0; // gene overexpressions,
  int del = 0;  // deletions
};

struct chip
{
  int chipNumber;                           // unique chip ID
  std::vector<sparse> perturbations;        // col2
  std::vector<std::string> perturbationsID; // col2
  int treatment;
  std::vector<std::string> deletedGenes;
  std::vector<std::string> overexpressedGenes;

  float time;
  std::vector<int> controls; // stores the indexes of the controls that pair
                             // with this measurement
  int totalRepeats = 1; // at least is always one
  std::vector<chip*> controlsPtr;
  int experiment; // col 1
  int replicaNumber;
};

bool
isPerturbationEqualSet(chip* c1, chip* c2);
bool
isPerturbed(int gene, chip* c);
bool
compareChips(chip c1, chip c2);

int
totalPerturbed(chip* c1, chip* c2);

int
readExpressionData(std::string filename);
void
readFeatures(std::string filename);
int
readGeneFile(std::string filename);
void
mapChipVector();
void
calculateMean();
void
selectPairs();
double
twoWayAnova(int genePair[2], float wko);
void
analizeGenePairs(std::string filename, float wko);
void
analizePartialPairs(std::vector<std::string> targets,
                    std::string filename,
                    float wko);

template<typename T>
int
countPerturbed(std::vector<T> v1, std::vector<T> v2)
{
  std::set<T> s1;
  std::set<T> s2;
  std::set<T> result;
  //  const bool isSparse = std::is_same<T, sparse>::value;
  for (T i : v1)
    s1.insert(i);
  for (T i : v2)
    s2.insert(i);
  std::set_difference(s1.begin(),
                      s1.end(),
                      s2.begin(),
                      s2.end(),
                      std::inserter(result, result.begin()));
  return result.size();
}
