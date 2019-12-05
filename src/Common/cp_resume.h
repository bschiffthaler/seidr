#pragma once

#include <sstream>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

template <typename T>
class cp_resume
{
public:
  cp_resume(T& param) {_param = param;}
  void print();
  bool save();
  bool resume();
private:
  T _param;
};

template <typename T>
void cp_resume<T>::print()
{
  std::stringstream ss;
  boost::archive::xml_oarchive oa(ss);
  oa << BOOST_SERIALIZATION_NVP(_param);
  std::cerr << ss.str();
}