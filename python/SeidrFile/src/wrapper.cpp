#include <boost/python.hpp>
#include <Serialize.h>
using namespace boost::python;

class SFWrapper
{
public:
  SFWrapper(const char * path);
  dict attr();
  list nodes();
  list algorithms();
  list supplementary();
  list betweenness();
  list closeness();
  list strength();
  list eigenvector();
  list pagerank();
  list katz();
  dict edge();
  void reset();
  bool next_edge();
  tuple edge_scores();
  tuple edge_ranks();
  dict edge_supplementary();
  tuple edge_index();
  uint8_t edge_flag();
private:
  SeidrFile _sf;
  SeidrFileHeader _sfh;
  SeidrFileEdge _e;
  uint32_t _i;
  uint32_t _j;
  uint64_t _k;
};

SFWrapper::SFWrapper(const char * path) :
  _sf(path), _i(0), _j(0), _k(0)
{
  _sf.open("r");
  _sfh.unserialize(_sf);
}

dict SFWrapper::attr()
{
  dict result;
  result["nodes"] = _sfh.attr.nodes;
  result["edges"] = _sfh.attr.edges;
  result["algorithms"] = _sfh.attr.nalgs;
  result["dense"] = _sfh.attr.dense;
  result["nsupp"] = _sfh.attr.nsupp;
  result["nsupp_str"] = _sfh.attr.nsupp_str;
  result["nsupp_int"] = _sfh.attr.nsupp_int;
  result["nsupp_flt"] = _sfh.attr.nsupp_flt;
  result["pagerank"] = _sfh.attr.pagerank_calc ? true : false;
  result["closeness"] = _sfh.attr.closeness_calc ? true : false;
  result["betweenness"] = _sfh.attr.betweenness_calc ? true : false;
  result["strength"] = _sfh.attr.strength_calc ? true : false;
  result["eigenvector"] = _sfh.attr.eigenvector_calc ? true : false;
  result["katz"] = _sfh.attr.katz_calc ? true : false;
  return result;
}

list SFWrapper::nodes()
{
  list result;
  for (const auto& n : _sfh.nodes)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::betweenness()
{
  list result;
  for (const auto& n : _sfh.betweenness)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::closeness()
{
  list result;
  for (const auto& n : _sfh.closeness)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::pagerank()
{
  list result;
  for (const auto& n : _sfh.pagerank)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::strength()
{
  list result;
  for (const auto& n : _sfh.strength)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::eigenvector()
{
  list result;
  for (const auto& n : _sfh.eigenvector)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::katz()
{
  list result;
  for (const auto& n : _sfh.katz)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::algorithms()
{
  list result;
  for (const auto& n : _sfh.algs)
  {
    result.append(n);
  }
  return result;
}

list SFWrapper::supplementary()
{
  list result;
  for (const auto& n : _sfh.supp)
  {
    result.append(n);
  }
  return result;
}

bool SFWrapper::next_edge()
{
  // Get next edge
  _e = SeidrFileEdge();
  while (_k < _sfh.attr.edges)
  {
    _e.unserialize(_sf, _sfh);
    if (_sfh.attr.dense)
    {
      _j++;
      if (_i <= _j)
      {
        _i++;
        _j = 0;
      }
      _e.index.i = _i;
      _e.index.j = _j;
    }
    else
    {
      _i = _e.index.i;
      _j = _e.index.j;
    }
    if (EDGE_EXISTS(_e.attr.flag))
    {
      _k++;
      return true;
    }
  }
  return false;
}

dict SFWrapper::edge()
{
  dict result;
  result.clear();
  if (_k <= _sfh.attr.edges && EDGE_EXISTS(_e.attr.flag))
  {
    if (! _sfh.attr.dense)
    {
      _i = _e.index.i;
      _j = _e.index.j;
    }
    result["flag"] = _e.attr.flag;
    result["index"] = make_tuple(_i, _j);
    list scores;
    list ranks;
    list supp_str;
    list supp_int;
    list supp_flt;
    for (const auto& s : _e.scores)
    {
      scores.append(s.s);
      ranks.append(s.r);
    }
    for (const auto& s : _e.supp_str)
    {
      supp_str.append(s);
    }
    for (const auto& s : _e.supp_int)
    {
      supp_int.append(s);
    }
    for (const auto& s : _e.supp_flt)
    {
      supp_flt.append(s);
    }
    result["scores"] = scores;
    result["ranks"] = ranks;
    result["supp_str"] = supp_str;
    result["supp_int"] = supp_int;
    result["supp_flt"] = supp_flt;
  }
  return result;
}

tuple SFWrapper::edge_scores()
{
  tuple result;
  if (_k <= _sfh.attr.edges && EDGE_EXISTS(_e.attr.flag))
  {
    list scores;
    for (const auto& s : _e.scores)
    {
      scores.append(s.s);
    }
    result = make_tuple(_e.index.i, _e.index.j, scores);
  }
  return result;
}

tuple SFWrapper::edge_ranks()
{
  tuple result;
  if (_k <= _sfh.attr.edges && EDGE_EXISTS(_e.attr.flag))
  {
    list scores;
    for (const auto& s : _e.scores)
    {
      scores.append(s.r);
    }
    result = make_tuple(_e.index.i, _e.index.j, scores);
  }
  return result;
}

tuple SFWrapper::edge_index()
{
  if (_k >  0 && _k <= _sfh.attr.edges)
  {
    return make_tuple(_e.index.i, _e.index.j);
  }
  return tuple();
}

uint8_t SFWrapper::edge_flag()
{
  uint8_t ret = 0;
  if (_k <= _sfh.attr.edges)
  {
    ret = _e.attr.flag;
  }
  return ret;
}

dict SFWrapper::edge_supplementary()
{
  dict result;
  result.clear();
  if (_k <= _sfh.attr.edges && EDGE_EXISTS(_e.attr.flag))
  {
    list supp_str;
    list supp_int;
    list supp_flt;
    for (const auto& s : _e.supp_str)
    {
      supp_str.append(s);
    }
    for (const auto& s : _e.supp_int)
    {
      supp_int.append(s);
    }
    for (const auto& s : _e.supp_flt)
    {
      supp_flt.append(s);
    }
    result["str"] = supp_str;
    result["int"] = supp_int;
    result["flt"] = supp_flt;
  }
  return result;
}


void SFWrapper::reset()
{
  _sf.seek(0);
  _sfh.unserialize(_sf);
  _i = 0;
  _j = 0;
  _k = 0;
}

BOOST_PYTHON_MODULE(libseidrfile)
{
  class_<SFWrapper>("SeidrFile_Internal", init<const char *>())
  .def("attr", &SFWrapper::attr)
  .def("nodes", &SFWrapper::nodes)
  .def("algorithms", &SFWrapper::algorithms)
  .def("supplementary", &SFWrapper::supplementary)
  .def("betweenness", &SFWrapper::betweenness)
  .def("closeness", &SFWrapper::closeness)
  .def("strength", &SFWrapper::strength)
  .def("eigenvector", &SFWrapper::eigenvector)
  .def("pagerank", &SFWrapper::pagerank)
  .def("katz", &SFWrapper::katz)
  .def("edge", &SFWrapper::edge)
  .def("reset", &SFWrapper::reset)
  .def("next_edge", &SFWrapper::next_edge)
  .def("edge_scores", &SFWrapper::edge_scores)
  .def("edge_ranks", &SFWrapper::edge_ranks)
  .def("edge_index", &SFWrapper::edge_index)
  .def("edge_supplementary", &SFWrapper::edge_supplementary)
  .def("edge_flag", &SFWrapper::edge_flag)
  ;
}
