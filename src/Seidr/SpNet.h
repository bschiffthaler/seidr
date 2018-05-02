#pragma once

#include <common.h>
#include <vector>

template <typename T>
class SpNet{
 public:
  size_t N;
  std::vector<T> map;
  std::vector<T> scores;
  std::vector<char> direct;
  SpNet(std::vector<size_t>&,
	std::vector<T>&,
	std::vector<char>&,
        std::vector<T>&,
	size_t);
  std::pair<T, char> get_by_coord(size_t, size_t);
  std::pair<T, char> get_by_index(size_t);
  T get_score_by_index(size_t);
  std::pair<size_t, size_t> in(size_t);
  size_t size();
};

// Constructor to insert values into the map array.
template <typename T>
SpNet<T>::SpNet(std::vector<size_t>& u,
                std::vector<T>& v,
                std::vector<char>& e,
                std::vector<T>& s,
                size_t n)
{
  N = n;
  /* We calculate the maximum number of values a lower triangular
     matrix needs to hold considering the number of genes. We
     also calculate the maximum rank for missing edges based on
     how many total edges we observe.
  */
  size_t lt_size = ((n * (n - 1)) / 2);
  T nan = std::numeric_limits<double>::signaling_NaN();
  map.resize(lt_size);
  map.assign(lt_size, nan);
  direct.resize(lt_size);
  direct.assign(lt_size, nan);
  scores.resize(lt_size);
  scores.assign(lt_size, nan);
  size_t vi = 0;
  for(size_t i = 0; i < u.size();){
    size_t row = u[i++] + 1;
    size_t col = u[i++] + 1;

    size_t diag = row * (row - 1) / 2;
    size_t diff = row - col;

    size_t index = diag - diff;
    map.at(index) = v[vi];
    direct.at(index) = e[vi];
    scores.at(index) = s[vi];
    vi++;
  }
}

/* Get a value by 2D coordiantes i and j.
   Time: O(1).*/
template <typename T>
std::pair<T, char> SpNet<T>::get_by_coord(size_t i, size_t j)
{
  size_t row = i + 1;
  size_t col = j + 1;

  size_t diag = row * (row - 1) / 2;
  size_t diff = row - col;

  size_t index = diag - diff;
  std::pair<T, char> ret(map[index], direct[index]);
  return ret;
}

/* Get a value by its array index. Time O(1).*/
template <typename T>
std::pair<T, char> SpNet<T>::get_by_index(size_t n)
{
  std::pair<T, char> ret(map[n], direct[n]);
  return ret;
}

template <typename T>
T SpNet<T>::get_score_by_index(size_t n)
{
  return scores[n];
}

/* Convert an array index to be a 2D coordiante in
 the lower triangular matrix. Time O(N) where N is
 the number of genes.*/
template <typename T>
std::pair<size_t, size_t> SpNet<T>::in(size_t n)
{
  size_t i = 1;
  size_t diag = 0;
  while(diag < n)
    {
      diag += ++i;
    }
  size_t diff = diag - n;
  size_t j = i - diff - 1;
  return std::pair<size_t, size_t>(i,j);
}

template <typename T>
size_t SpNet<T>::size(){
  return N;
}

