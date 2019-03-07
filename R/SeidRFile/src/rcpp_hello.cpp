#include <Rcpp.h>
#include <cmath>
#include "../inst/include/SeidRFile.h"
using namespace Rcpp;

// [[Rcpp::export]]
XPtr<SeidrFile> SeidrFile__ptr(std::string p) {
  SeidrFile * ptr = new SeidrFile(p.c_str());
  ptr->open("r");
  return XPtr<SeidrFile>(ptr);
}

// [[Rcpp::export]]
XPtr<SeidrFileHeader> SeidrFileHeader__ptr(SEXP xp) {
  XPtr<SeidrFile> ptr(xp);
  SeidrFileHeader * h = new SeidrFileHeader();
  h->unserialize(*ptr);
  return XPtr<SeidrFileHeader>(h);
}

// [[Rcpp::export]]
List SeidrFileHeader__ptr__getAttr(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  uint64_t nodes = ptr->attr.nodes;
  uint64_t edges = ptr->attr.edges;
  uint16_t nalgs = ptr->attr.nalgs;
  uint8_t  dense = ptr->attr.dense;
  uint32_t nsupp = ptr->attr.nsupp;
  uint16_t nsupp_str = ptr->attr.nsupp_str;
  uint16_t nsupp_int = ptr->attr.nsupp_int;
  uint16_t nsupp_flt = ptr->attr.nsupp_flt;
  uint8_t  pagerank_calc = ptr->attr.pagerank_calc;
  uint8_t  closeness_calc = ptr->attr.closeness_calc;
  uint8_t  betweenness_calc = ptr->attr.betweenness_calc;
  uint8_t  strength_calc = ptr->attr.strength_calc;
  uint8_t  eigenvector_calc = ptr->attr.eigenvector_calc;
  uint8_t  katz_calc = ptr->attr.katz_calc;
  List attr = List::create(nodes, edges, nalgs, dense, nsupp,
                           nsupp_str, nsupp_int, nsupp_flt,
                           pagerank_calc, closeness_calc,
                           betweenness_calc, strength_calc,
                           eigenvector_calc, katz_calc);
  return attr;
}

// [[Rcpp::export]]
SEXP SeidrFileHeader__ptr__Position(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return List::create(ptr->ci, ptr->cj, ptr->k);
}

// [[Rcpp::export]]
CharacterVector SeidrFileHeader__ptr__getNodes(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->nodes);
}

// [[Rcpp::export]]
CharacterVector SeidrFileHeader__ptr__getAlgs(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->algs);
}

// [[Rcpp::export]]
CharacterVector SeidrFileHeader__ptr__getSTags(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->supp);
}

// [[Rcpp::export]]
NumericVector SeidrFileHeader__ptr__getPageRank(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->pagerank);
}

// [[Rcpp::export]]
NumericVector SeidrFileHeader__ptr__getCloseness(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->closeness);
}

// [[Rcpp::export]]
NumericVector SeidrFileHeader__ptr__getBetweenness(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->betweenness);
}

// [[Rcpp::export]]
NumericVector SeidrFileHeader__ptr__getStrength(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->strength);
}

// [[Rcpp::export]]
NumericVector SeidrFileHeader__ptr__getEigenvector(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->eigenvector);
}

// [[Rcpp::export]]
NumericVector SeidrFileHeader__ptr__getKatz(SEXP xp) {
  XPtr<SeidrFileHeader> ptr(xp);
  return wrap(ptr->katz);
}

// [[Rcpp::export]]
SEXP SeidrFile__ptr__vectorizeSF(SEXP xpf, SEXP xph, bool scores, bool ranks, bool flags,
                                 bool supp_str, bool supp_int, bool supp_flt,
                                 bool indices, uint64_t chunksize) {
  // Get smart pointers to data structures
  XPtr<SeidrFile> fptr(xpf);
  XPtr<SeidrFileHeader> hptr(xph);

  // If we are not chunking get all
  if (chunksize == 0) chunksize = hptr->attr.edges;

  // R return data destinations
  std::vector<std::vector<double>> rets;
  std::vector<std::vector<double>> retr;
  std::vector<uint32_t> is;
  std::vector<uint32_t> js;
  std::vector<uint8_t> fs;
  std::vector<std::vector<std::string>> sups;
  std::vector<std::vector<int>> supi;
  std::vector<std::vector<float>> supf;

  // Header index data
  uint32_t& i = hptr->ci;
  uint32_t& j = hptr->cj;
  uint64_t& k = hptr->k;
  uint64_t read = 0;

  if (hptr->attr.dense) {
    while(k < hptr->attr.edges && read < chunksize) {
      SeidrFileEdge e;
      e.unserialize((*fptr), (*hptr));
      if (i == 0) { // First edge
        i = 1; j = 0;
      }
      else if (j == (i - 1)) { // subsequent edges
        i++; j = 0;
      } else {
        j++;
      }
      if (EDGE_EXISTS(e.attr.flag)) {
        // Update index data
        e.index.i = i;
        e.index.j = j;
        hptr->ci = i;
        hptr->cj = j;
        hptr->k++;
        read++;
        // Get all user requested fields
        if (indices) {
          is.push_back(i + 1);
          js.push_back(j + 1);
        }
        if (flags) {
          fs.push_back(e.attr.flag);
        }

        if (supp_str) {
          std::vector<std::string> sup_s;
          for (auto& x : e.supp_str)
            sup_s.push_back(x);
          sups.push_back(sup_s);
        }

        if (supp_int) {
          std::vector<int> sup_i;
          for (auto& x : e.supp_int)
            sup_i.push_back(x);
          supi.push_back(sup_i);
        }

        if (supp_flt) {
          std::vector<float> sup_f;
          for (auto& x : e.supp_flt)
            sup_f.push_back(x);
          supf.push_back(sup_f);
        }

        std::vector<double> svs;
        std::vector<double> svr;
        for (auto& x : e.scores) {
          if (ranks) {
            if (std::isfinite(x.r)) {
              svr.push_back(x.r);
            } else {
              svr.push_back(NA_REAL);
            }
          }
          if (scores) {
            if (std::isfinite(x.s)) {
              svs.push_back(x.s);
            } else {
              svs.push_back(NA_REAL);
            }
          }
        }
        if (ranks) {
          retr.push_back(svr);
        }
        if (scores) {
          rets.push_back(svs);
        }
      }
    }
  } else {
    while(k < hptr->attr.edges && read < chunksize) {
      SeidrFileEdge e;
      e.unserialize((*fptr), (*hptr));

      // Update indices
      hptr->ci = e.index.i;
      hptr->cj = e.index.j;
      hptr->k++;
      read++;

      // Get fields user requested
      if (indices) {
        is.push_back(e.index.i + 1);
        js.push_back(e.index.j + 1);
      }

      if (flags) {
        fs.push_back(e.attr.flag);
      }

      if (supp_str) {
        std::vector<std::string> sup_s;
        for (auto& x : e.supp_str)
          sup_s.push_back(x);
        sups.push_back(sup_s);
      }

      if (supp_int) {
        std::vector<int> sup_i;
        for (auto& x : e.supp_int)
          sup_i.push_back(x);
        supi.push_back(sup_i);
      }

      if (supp_flt) {
        std::vector<float> sup_f;
        for (auto& x : e.supp_flt)
          sup_f.push_back(x);
        supf.push_back(sup_f);
      }

      std::vector<double> svs;
      std::vector<double> svr;
      for (auto& x : e.scores) {
        if (ranks) {
          if (std::isfinite(x.r)) {
            svr.push_back(x.r);
          } else {
            svr.push_back(NA_REAL);
          }
        }
        if (scores) {
          if (std::isfinite(x.s)) {
            svs.push_back(x.s);
          } else {
            svs.push_back(NA_REAL);
          }
        }
      }
      if (ranks) {
        retr.push_back(svr);
      }
      if (scores) {
        rets.push_back(svs);
      }
    }
  }
  return List::create(wrap(rets), wrap(retr), wrap(is), wrap(js), wrap(fs),
                      wrap(sups), wrap(supi), wrap(supf));
}

// [[Rcpp::export]]
SEXP SeidrFile__ptr__vectorizeSingle(SEXP xpf, SEXP xph, bool score, bool rank,
                                     bool indices, uint32_t index,
                                     uint64_t chunksize) {
  // Get smart pointers to data structures
  XPtr<SeidrFile> fptr(xpf);
  XPtr<SeidrFileHeader> hptr(xph);

  // 0 offset in C++
  index--;

  // If we are not chunking get all
  if (chunksize == 0) chunksize = hptr->attr.edges;

  // R return data destinations
  std::vector<double> rets;
  std::vector<double> retr;
  std::vector<uint32_t> is;
  std::vector<uint32_t> js;

  // Header index data
  uint32_t& i = hptr->ci;
  uint32_t& j = hptr->cj;
  uint64_t& k = hptr->k;
  uint64_t read = 0;

  if (hptr->attr.dense) {
    while(k < hptr->attr.edges && read < chunksize) {
      SeidrFileEdge e;
      e.unserialize((*fptr), (*hptr));
      if (i == 0) { // First edge
        i = 1; j = 0;
      }
      else if (j == (i - 1)) { // subsequent edges
        i++; j = 0;
      } else {
        j++;
      }
      if (EDGE_EXISTS(e.attr.flag)) {
        // Update index data
        e.index.i = i;
        e.index.j = j;
        hptr->ci = i;
        hptr->cj = j;
        hptr->k++;
        read++;
        // Get all user requested fields
        if (indices) {
          is.push_back(i + 1);
          js.push_back(j + 1);
        }

        if (score) {
          if (std::isfinite(e.scores[index].s)) {
            rets.push_back(e.scores[index].s);
          } else {
            rets.push_back(NA_REAL);
          }
        }

        if (rank) {
          if (std::isfinite(e.scores[index].r)) {
            retr.push_back(e.scores[index].r);
          } else {
            retr.push_back(NA_REAL);
          }
        }

      }
    }
  } else {
    while(k < hptr->attr.edges && read < chunksize) {
      SeidrFileEdge e;
      e.unserialize((*fptr), (*hptr));

      // Update indices
      hptr->ci = e.index.i;
      hptr->cj = e.index.j;
      hptr->k++;
      read++;

      // Get fields user requested
      if (indices) {
        is.push_back(e.index.i + 1);
        js.push_back(e.index.j + 1);
      }

      if (score) {
        if (std::isfinite(e.scores[index].s)) {
          rets.push_back(e.scores[index].s);
        } else {
          rets.push_back(NA_REAL);
        }
      }

      if (rank) {
        if (std::isfinite(e.scores[index].r)) {
          retr.push_back(e.scores[index].r);
        } else {
          retr.push_back(NA_REAL);
        }
      }
    }
  }
  return List::create(wrap(rets), wrap(retr), wrap(is), wrap(js));
}

// [[Rcpp::export]]
void SeidrFile__ptr__close(SEXP xpf) {
  XPtr<SeidrFile>(xpf)->close();
}

// [[Rcpp::export]]
SEXP SeidrFile__ptr__isopen(SEXP xpf) {
  return wrap(XPtr<SeidrFile>(xpf)->_opened);
}

