#pragma once

#include <vector>
#include <string>
#include <armadillo>
#include <map>
#include <fstream>

#ifdef SEIDR_SCORE_DOUBLE
typedef arma::mat mat_t;
#else
typedef arma::fmat mat_t;
#endif

void parse_el(mat_t& m,
              std::map<std::string, size_t>& gm,
              std::istream& ifs,
              char sep);
void parse_sm(mat_t& m, std::istream& ifs,
              char sep);
void parse_ltri(mat_t& m, std::istream& ifs, bool diag, char sep);
void parse_utri(mat_t& m, std::istream& ifs, bool diag, char sep);
void parse_aracne(mat_t& m, std::map<std::string, size_t>& mp,
                  std::istream& ifs, char sep);
void write_el(mat_t& m, std::vector<std::string>& g,
              seidr_score_t fill, char sep);
void write_sm(mat_t& m, char sep);
void write_ltri(mat_t& m, bool diag, char sep);
void write_utri(mat_t& m, bool diag, char sep);
void write_aracne(mat_t& m, std::vector<std::string>& g,
                  char sep);
int convert(int argc, char ** argv);
