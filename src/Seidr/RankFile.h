#pragma once

#include <string>
#include <vector>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

class RankFileHeader {
private:
  size_t ne;
  size_t nn;
  char type;
  std::string v;
  std::string cmd;
  std::vector<std::string> node_names;
  std::vector<std::string> method_names;
  char sort_order;
public:
  friend class cereal::access;
  template <class Archive>
  void serialize( Archive & ar )
    {
      ar(ne, nn, v, cmd, node_names, type,
	 method_names, sort_order);
    }
  void set_ne(size_t x);
  void set_nn(size_t x);
  void set_v(std::string x);
  void set_cmd(std::string x);
  void set_node_names(std::vector<std::string> x);
  void set_method_names(std::vector<std::string> x);
  void set_type(int x);
  void set_sort_order(char);
  size_t get_ne();
  size_t get_nn();
  std::string get_v();
  std::string get_cmd();
  std::vector<std::string> get_node_names();
  std::vector<std::string> get_method_names();
  char get_type();
  char get_sort_order();
};

template <class T>
class RankFileEdge {
private:
  size_t i; //from
  size_t j; //to
  T r; //rank
  std::vector<T> w; //weight
  char d; //directionality
public:
  friend class cereal::access;
  template <class Archive>
  void serialize( Archive & ar )
    {
      ar(i, j, w, r, d);
    }
  void set_i(size_t);
  void set_j(size_t);
  void set_w(T, size_t);
  void push_back_w(T);
  void set_r(T);
  void set_d(char);
  size_t get_i();
  size_t get_j();
  std::vector<T> get_w();
  T get_w_at(size_t);
  T get_r();
  char get_d();
};

inline void RankFileHeader::set_ne(size_t x) {ne = x;}
inline void RankFileHeader::set_nn(size_t x) {nn = x;}
inline void RankFileHeader::set_v(std::string x) {v = x;}
inline void RankFileHeader::set_cmd(std::string x) {cmd = x;}
inline void RankFileHeader::set_node_names(std::vector<std::string> x) {node_names = x;}
inline void RankFileHeader::set_method_names(std::vector<std::string> x) {method_names = x;}
inline void RankFileHeader::set_type(int x) {type = x;}
inline size_t RankFileHeader::get_ne(){return ne;}
inline size_t RankFileHeader::get_nn(){return nn;}
inline std::string RankFileHeader::get_v(){return v;}
inline std::string RankFileHeader::get_cmd(){return cmd;}
inline std::vector<std::string> RankFileHeader::get_node_names(){return node_names;}
inline std::vector<std::string> RankFileHeader::get_method_names(){return method_names;}
inline char RankFileHeader::get_type(){return type;}
inline char RankFileHeader::get_sort_order(){return sort_order;}

template <class T>
  void RankFileEdge<T>::set_i(size_t x) {i = x;}

template <class T>
  void RankFileEdge<T>::set_j(size_t x) {j =x;}

template <class T>
  void RankFileEdge<T>::set_w(T x, size_t id) {w[id] = x;}

template <class T>
void RankFileEdge<T>::push_back_w(T x) {w.push_back(x);}

template <class T>
  void RankFileEdge<T>::set_r(T x) {r = x;}

template <class T>
  void RankFileEdge<T>::set_d(char x) {d = x;}

template <class T>
  size_t RankFileEdge<T>::get_i() {return i;}

template <class T>
  size_t RankFileEdge<T>::get_j() {return j;}

template <class T>
std::vector<T> RankFileEdge<T>::get_w() {return w;}

template <class T>
T RankFileEdge<T>::get_w_at(size_t pos) {return w[pos];}

template <class T>
  T RankFileEdge<T>::get_r() {return r;}

template <class T>
  char RankFileEdge<T>::get_d() {return d;}

