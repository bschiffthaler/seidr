#pragma once

#include <string>
#include <vector>
#include <set>
#include <armadillo>
#include <algorithm>

typedef arma::uword seidr_index;

struct sparse {

  std::string index; //perturbation ID
  float value=0; //weight of the perturbation

  sparse(std::string idx, float val): index(idx), value(val){}

  bool operator<(const sparse& s2) const {
    if (index < s2.index) {
      return true;
    } else if (index == s2.index) {
      if(value < s2.value)
	return true;
    }
    return false;
  }
  bool operator==(const sparse &s2) const {
    return (index == s2.index) && (value == s2.value);
  }
};

struct gene {
  int over=0; //gene overexpressions, 
  int del=0; // deletions
};

struct chip {
  int chipNumber; //unique chip ID
  std::vector<sparse> perturbations; //col2
  std::vector<std::string> perturbationsID; //col2
  int treatment;
  std::vector<std::string> deletedGenes;
  std::vector<std::string> overexpressedGenes;
  
  float time;
  std::vector<int> controls; //stores the indexes of the controls that pair with this measurement
  int totalRepeats=1;// at least is always one
  std::vector<chip*> controlsPtr;
  int experiment; //col 1
  int replicaNumber;
};


bool isPerturbationEqualSet (chip* c1, chip* c2);
bool isPerturbed(int gene, chip* c);
bool compareChips(chip c1, chip c2);

int totalPerturbed(chip* c1, chip* c2);

void readExpressionData(std::string filename);
void readFeatures(std::string filename);
void readGeneFile(std::string filename);
void mapChipVector();
void calculateMean();
void selectPairs();
double twoWayAnova (int genePair[2], float wko);
void analizeGenePairs(std::string filename, float wko);
void analizePartialPairs(std::vector<std::string> targets,
			 std::string filename, float wko);

template<typename T>
int countPerturbed(std::vector<T> v1, std::vector<T> v2) {
  std::set<T> s1;
  std::set<T> s2;
  std::set<T> result;
  //  const bool isSparse = std::is_same<T, sparse>::value;
  for (T i:v1)
    s1.insert(i);
  for (T i:v2)
    s2.insert(i);
  std::set_difference (s1.begin(), s1.end(), s2.begin(), s2.end(),
                       std::inserter(result, result.begin()) );
  return result.size();
}
