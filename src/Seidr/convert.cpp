//
// Seidr - Create and operate on gene crowd networks
// Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
//
// This file is part of Seidr.
//
// Seidr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Seidr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
//

#include <BSlogger.hpp>
#include <common.h>
#include <convert.h>
#include <fs.h>

#include <armadillo>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace po = boost::program_options;

void
parse_el(mat_t& m,
         std::map<std::string, size_t>& gm,
         std::istream& ifs,
         char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  while (std::getline(ifs, line)) {
    std::string field;
    std::string from;
    std::string to;
    seidr_score_t w = 0;
    unsigned int counter = 0;
    std::stringstream ss(line);
    while (std::getline(ss, field, sep)) {
      if (counter == 0) {
        from = field;
      }
      if (counter == 1) {
        to = field;
      }
      if (counter == 2) {
        w = std::stod(field);
      }
      counter++;
    }
    size_t i = gm[from];
    size_t j = gm[to];
    m(i, j) = w;
  }
}

void
parse_sm(mat_t& m, std::istream& ifs, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  size_t i = 0;
  while (std::getline(ifs, line)) {
    std::string field;
    std::stringstream ss(line);
    size_t j = 0;
    while (std::getline(ss, field, sep)) {
      seidr_score_t w = 0;
      w = std::stod(field);
      m(i, j) = w;
      j++;
    }
    i++;
  }
}

void
parse_ltri(mat_t& m, std::istream& ifs, bool diag, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  size_t i = diag ? 0 : 1;
  while (std::getline(ifs, line)) {
    size_t j = 0;
    std::string field;
    seidr_score_t w = 0;
    std::stringstream ss(line);
    while (std::getline(ss, field, sep)) {
      w = std::stod(field);
      m(i, j) = w;
      m(j, i) = w;
      j++;
    }
    i++;
  }
}

void
parse_utri(mat_t& m, std::istream& ifs, bool diag, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  size_t i = 0;
  while (std::getline(ifs, line)) {
    size_t j = diag ? i : (i + 1);
    std::string field;
    seidr_score_t w = 0;
    std::stringstream ss(line);
    ss.ignore(i);
    while (std::getline(ss, field, sep)) {
      w = std::stod(field);
      m(i, j) = w;
      m(j, i) = w;
      j++;
    }
    i++;
  }
}

void
parse_aracne(mat_t& m,
             std::map<std::string, size_t>& mp,
             std::istream& ifs,
             char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  while (std::getline(ifs, line)) {
    if (line[0] == '>') {
      continue;
    }

    std::string field;
    seidr_score_t w = 0;
    size_t i = 0;
    size_t j = 0;

    std::stringstream ss(line);
    std::getline(ss, field, sep);
    i = mp[field];

    while (std::getline(ss, field, sep)) {
      j = mp[field];
      std::getline(ss, field, sep);
      w = std::stod(field);
      m(i, j) = w;
    }
  }
}

void
write_el(mat_t& m, std::vector<std::string>& g, std::ostream& out, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  for (size_t i = 0; i < g.size(); i++) {
    for (size_t j = 0; j < g.size(); j++) {
      if (arma::is_finite(m(i, j))) {
        out << g[i] << sep << g[j] << sep << m(i, j) << '\n';
      }
    }
  }
}

void
write_sm(mat_t& m, std::ostream& out, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  for (size_t i = 0; i < m.n_rows; i++) {
    for (size_t j = 0; j < m.n_cols; j++) {
      out << m(i, j) << (j == m.n_cols - 1 ? '\n' : sep);
    }
  }
}

void
write_ltri(mat_t& m, bool diag, std::ostream& out, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  if (diag) {
    for (size_t i = 0; i < m.n_rows; i++) {
      for (size_t j = 0; j <= i; j++) {
        out << m(i, j) << (j == i ? '\n' : sep);
      }
    }
  } else {
    for (size_t i = 1; i < m.n_rows; i++) {
      for (size_t j = 0; j < i; j++) {
        out << m(i, j) << (j == i - 1 ? '\n' : sep);
      }
    }
  }
}

void
write_utri(mat_t& m, bool diag, std::ostream& out, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  if (diag) {
    for (size_t i = 0; i < m.n_rows; i++) {
      for (size_t j = 0; j < i; j++) {
        out << sep;
      }
      for (size_t j = i; j < m.n_cols; j++) {
        out << m(i, j) << (j == m.n_cols - 1 ? '\n' : sep);
      }
    }
  } else {
    for (size_t i = 0; i < m.n_rows - 1; i++) {
      for (size_t j = 0; j < i; j++) {
        out << sep;
      }
      for (size_t j = i + 1; j < m.n_cols; j++) {
        out << m(i, j) << (j == m.n_cols - 1 ? '\n' : sep);
      }
    }
  }
}

void
write_aracne(mat_t& m, std::vector<std::string>& g, std::ostream& out, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  for (size_t i = 0; i < m.n_rows; i++) {
    out << g[i] << sep;
    for (size_t j = 0; j < m.n_cols; j++) {
      out << g[j] << sep;
      out << m(i, j) << (j == m.n_cols - 1 ? '\n' : sep);
    }
  }
}

int
convert(const std::vector<std::string>& args)
{

  logger log(std::cerr, "convert");

  try {

    seidr_conv_param_t param;

    po::options_description umbrella("Convert different text based formats");

    po::options_description opt("Common Options");
    opt.add_options()("force,f",
                      po::bool_switch(&param.force)->default_value(false),
                      "Force overwrite output file if it exists")(
      "help,h", "Show this help message");

    po::options_description iopt("Input Options");
    iopt.add_options()(
      "infile,i",
      po::value<std::string>(&param.infile)->default_value("-"),
      "Input file name ['-' for stdin]")(
      "genes,g",
      po::value<std::string>(&param.gene_file)->required(),
      "Input gene file name")(
      "from,q",
      po::value<std::string>(&param.in_format)->required(),
      "Input file format [edge-list, sym-mat, low-tri, up-tri, low-tri-diag, "
      "up-tri-diag, aracne]")(
      "in-separator,s",
      po::value<std::string>(&param.in_sep)->default_value("", "\\t"),
      "Input separator");

    po::options_description oopt("Output Options");
    oopt.add_options()(
      "outfile,o",
      po::value<std::string>(&param.outfile)->default_value("-"),
      "Output file name ['-' for stdout]")(
      "to,t",
      po::value<std::string>(&param.out_format)->required(),
      "Output file format [edge-list, sym-mat, low-tri, up-tri, low-tri-diag, "
      "up-tri-diag, aracne]")(
      "out-separator,S",
      po::value<std::string>(&param.out_sep)->default_value("", "\\t"),
      "Output separator");

    po::options_description fopt("Format Options");
    fopt.add_options()(
      "fill,F",
      po::value<std::string>(&param.fill)->default_value("NaN"),
      "Fill value for missing data (Number, NaN, Inf, -Inf)")(
      "precision,p",
      po::value<uint16_t>(&param.prec)->default_value(SEIDR_CONVERT_DEF_PREC),
      "Number of decimals to print");

    umbrella.add(iopt).add(oopt).add(fopt).add(opt);

    po::variables_map vm;
    po::store(po::command_line_parser(args).options(umbrella).run(), vm);

    if (vm.count("help") != 0 || args.empty()) {
      std::cerr << umbrella << '\n';
      return 1;
    }

    po::notify(vm);

    if (param.infile != "-") {
      param.infile = to_absolute(param.infile);
      assert_exists(param.infile);
      assert_can_read(param.infile);
      assert_no_cr(param.infile);
    }

    if (param.outfile != "-") {
      param.outfile = to_absolute(param.outfile);
      if (!param.force) {
        assert_no_overwrite(param.outfile);
      }
      assert_dir_is_writeable(dirname(param.outfile));
    }
    assert_arg_constraint<std::string>({ "edge-list",
                                         "sym-mat",
                                         "low-tri",
                                         "up-tri",
                                         "low-tri-diag",
                                         "up-tri-diag",
                                         "aracne" },
                                       param.in_format);
    assert_arg_constraint<std::string>({ "edge-list",
                                         "sym-mat",
                                         "low-tri",
                                         "up-tri",
                                         "low-tri-diag",
                                         "up-tri-diag",
                                         "aracne" },
                                       param.out_format);

    if (vm["in-separator"].defaulted()) {
      param.in_sep = "\t";
    }

    if (vm["out-separator"].defaulted()) {
      param.out_sep = "\t";
    }

    std::vector<std::string> genes = read_genes(param.gene_file, '\n', '\t');
    seidr_score_t fi = 0;

    if (param.fill == "Inf") {
      fi = arma::datum::inf;
    } else if (param.fill == "-Inf") {
      fi = -arma::datum::inf;
    } else if (param.fill == "NaN") {
      fi = arma::datum::nan;
    } else {
      fi = std::stod(param.fill);
    }

    mat_t gm(genes.size(), genes.size());

    gm.fill(fi);

    std::map<std::string, size_t> gene_map;
    for (size_t i = 0; i < genes.size(); i++) {
      gene_map[genes[i]] = i;
    }

    std::shared_ptr<std::istream> in_stream;
    if (param.infile == "-") {
      in_stream.reset(&(std::cin), no_delete);
    } else {
      in_stream =
        std::shared_ptr<std::istream>(new std::ifstream(param.infile.c_str()));
    }

    std::shared_ptr<std::ostream> out_stream;
    if (param.outfile == "-") {
      out_stream.reset(&std::cout, no_delete);
    } else {
      out_stream =
        std::shared_ptr<std::ostream>(new std::ofstream(param.outfile.c_str()));
    }

    out_stream->precision(param.prec);
    out_stream->setf(std::ios::fixed, std::ios::floatfield);

    // Generic reading/writing if no specialized function is available
    if (param.in_format == "edge-list") {
      parse_el(gm, gene_map, *in_stream, param.in_sep[0]);
    } else if (param.in_format == "sym-mat") {
      parse_sm(gm, *in_stream, param.in_sep[0]);
    } else if (param.in_format == "low-tri") {
      parse_ltri(gm, *in_stream, false, param.in_sep[0]);
    } else if (param.in_format == "up-tri") {
      parse_utri(gm, *in_stream, false, param.in_sep[0]);
    } else if (param.in_format == "low-tri-diag") {
      parse_ltri(gm, *in_stream, true, param.in_sep[0]);
    } else if (param.in_format == "up-tri-diag") {
      parse_utri(gm, *in_stream, true, param.in_sep[0]);
    } else if (param.in_format == "aracne") {
      parse_aracne(gm, gene_map, *in_stream, param.in_sep[0]);
    } else {
      throw std::runtime_error("Unknown input format");
    }

    // Output
    if (param.out_format == "edge-list") {
      write_el(gm, genes, *out_stream, param.out_sep[0]);
    } else if (param.out_format == "sym-mat") {
      write_sm(gm, *out_stream, param.out_sep[0]);
    } else if (param.out_format == "low-tri") {
      write_ltri(gm, false, *out_stream, param.out_sep[0]);
    } else if (param.out_format == "up-tri") {
      write_utri(gm, false, *out_stream, param.out_sep[0]);
    } else if (param.out_format == "low-tri-diag") {
      write_ltri(gm, true, *out_stream, param.out_sep[0]);
    } else if (param.out_format == "up-tri-diag") {
      write_utri(gm, true, *out_stream, param.out_sep[0]);
    } else if (param.out_format == "aracne") {
      write_aracne(gm, genes, *out_stream, param.out_sep[0]);
    } else {
      throw std::runtime_error("Unknown output format");
    }

  } catch (const po::error& e) {
    log(LOG_ERR) << "[Argument Error]: " << e.what() << '\n';
    return 1;
  } catch (const std::runtime_error& e) {
    log(LOG_ERR) << "[Runtime Error]: " << e.what() << '\n';
    return 1;
  } catch (const std::exception& e) {
    log(LOG_ERR) << "[Generic Error]: " << e.what() << '\n';
    return 1;
  }

  return 0;
}
