#pragma once

#include <memory>
#include <cstdint>

template <typename T>
class LowerTriangular{
 public:
  LowerTriangular();
  LowerTriangular(uint64_t size);
  T&  at(uint64_t i, uint64_t j);
  const T& at(uint64_t i, uint64_t j) const;
  T get(uint64_t i, uint64_t j);
  const T get(uint64_t i, uint64_t j) const;
  uint64_t size(void) const {return _size;}
  void reset(uint64_t size);
 private:
  std::unique_ptr<T[]> _data;
  uint64_t _size;
};

template <typename T>
LowerTriangular<T>::LowerTriangular() :
_data(nullptr), _size(0)
{}

template <typename T>
LowerTriangular<T>::LowerTriangular(uint64_t size) :
_size(size)
{
  _data.reset(new T[size * (size - 1) / 2]);
}

template <typename T>
void LowerTriangular<T>::reset(uint64_t size)
{
  _size = size;
  _data.reset(new T[size * (size - 1) / 2]);
}

template<typename T>
T& LowerTriangular<T>::at(uint64_t i, uint64_t j)
{
  uint64_t row = i + 1;
  uint64_t col = j + 1;
  
  uint64_t diag = row * (row - 1) / 2;
  uint64_t diff = row - col;

  uint64_t index = diag - diff;
  return _data[index];
}

template<typename T>
const T& LowerTriangular<T>::at(uint64_t i, uint64_t j) const
{
  uint64_t row = i + 1;
  uint64_t col = j + 1;
  
  uint64_t diag = row * (row - 1) / 2;
  uint64_t diff = row - col;

  uint64_t index = diag - diff;
  const T& ret =  _data[index];
  return ret;
}

template<typename T>
T LowerTriangular<T>::get(uint64_t i, uint64_t j)
{
  uint64_t row = i + 1;
  uint64_t col = j + 1;

  uint64_t diag = row * (row - 1) / 2;
  uint64_t diff = row - col;

  uint64_t index = diag - diff;
  T ret =  (* _data[index]);
  return ret;
}

template<typename T>
const T LowerTriangular<T>::get(uint64_t i, uint64_t j) const
{
  uint64_t row = i + 1;
  uint64_t col = j + 1;

  uint64_t diag = row * (row - 1) / 2;
  uint64_t diff = row - col;

  uint64_t index = diag - diff;
  const T ret =  (* _data[index]);
  return ret;
}
