#include <BSlogger.h>
#include <common.h>
#include <fs.h>
#include <convert.h>
#include <tclap_format.h>

#include <iostream>
#include <vector>
#include <string>
#include <armadillo>
#include <map>
#include <fstream>
#include <sstream>
#include <tclap/CmdLine.h>
#include <memory>

void parse_el(mat_t& m,
              std::map<std::string, size_t>& gm,
              std::istream& ifs,
              char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  while (std::getline(ifs, line))
  {
    std::string field, from, to;
    seidr_score_t w = 0;
    unsigned int counter = 0;
    std::stringstream ss(line);
    while (std::getline(ss, field, sep))
    {
      if (counter == 0) from = field;
      if (counter == 1) to = field;
#ifdef SEIDR_SCORE_DOUBLE
      if (counter == 2) w = std::stod(field);
#else
      if (counter == 2) w = std::stof(field);
#endif
      counter++;
    }
    size_t i = gm[from];
    size_t j = gm[to];
    m(i, j) = w;
  }
}

void parse_sm(mat_t& m, std::istream& ifs,
              char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  size_t i = 0;
  while (std::getline(ifs, line))
  {
    std::string field;
    std::stringstream ss(line);
    size_t j = 0;
    while (std::getline(ss, field, sep))
    {
      seidr_score_t w;
#ifdef SEIDR_SCORE_DOUBLE
      w = std::stod(field);
#else
      w = std::stof(field);
#endif
      m(i, j) = w;
      j++;
    }
    i++;
  }
}

void parse_ltri(mat_t& m, std::istream& ifs, bool diag, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  size_t i = diag ? 0 : 1;
  while (std::getline(ifs, line))
  {
    size_t j = 0;
    std::string field;
    seidr_score_t w;
    std::stringstream ss(line);
    while (std::getline(ss, field, sep))
    {
#ifdef SEIDR_SCORE_DOUBLE
      w = std::stod(field);
#else
      w = std::stof(field);
#endif
      m(i, j) = w;
      m(j, i) = w;
      j++;
    }
    i++;
  }
}

void parse_utri(mat_t& m, std::istream& ifs, bool diag, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  size_t i = 0;
  while (std::getline(ifs, line))
  {
    size_t j = diag ? i : (i + 1);
    std::string field;
    seidr_score_t w;
    std::stringstream ss(line);
    ss.ignore(i);
    while (std::getline(ss, field, sep))
    {
#ifdef SEIDR_SCORE_DOUBLE
      w = std::stod(field);
#else
      w = std::stof(field);
#endif
      m(i, j) = w;
      m(j, i) = w;
      j++;
    }
    i++;
  }
}

void parse_aracne(mat_t& m, std::map<std::string, size_t>& mp,
                  std::istream& ifs, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  std::string line;
  while (std::getline(ifs, line))
  {
    if (line[0] == '>') continue;

    std::string field;
    seidr_score_t w;
    size_t i, j;

    std::stringstream ss(line);
    std::getline(ss, field, sep);
    i = mp[field];

    while (std::getline(ss, field, sep))
    {
      j = mp[field];
      std::getline(ss, field, sep);
#ifdef SEIDR_SCORE_DOUBLE
      w = std::stod(field);
#else
      w = std::stof(field);
#endif
      m(i, j) = w;
    }
  }
}

void write_el(mat_t& m, std::vector<std::string>& g,
              seidr_score_t fill, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  for (size_t i = 0; i < g.size(); i++)
  {
    for (size_t j = 0; j < g.size(); j++)
    {
      if (arma::is_finite(m(i, j)))
      {
        std::cout << g[i] << sep
                  << g[j] << sep
                  << m(i, j) << '\n';
      }
    }
  }
}

void write_sm(mat_t& m, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  for (size_t i = 0; i < m.n_rows; i++)
  {
    for (size_t j = 0; j < m.n_cols; j++)
    {
      std::cout << m(i, j)
                << ( j == m.n_cols - 1 ? '\n' : sep);
    }
  }
}

void write_ltri(mat_t& m, bool diag, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  if (diag)
  {
    for (size_t i = 0; i < m.n_rows; i++)
    {
      for (size_t j = 0; j <= i; j++)
      {
        std::cout << m(i, j)
                  << ( j == i ? '\n' : sep);
      }
    }
  }
  else
  {
    for (size_t i = 1; i < m.n_rows; i++)
    {
      for (size_t j = 0; j < i; j++)
      {
        std::cout << m(i, j)
                  << ( j == i - 1 ? '\n' : sep);
      }
    }
  }
}

void write_utri(mat_t& m, bool diag, char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  if (diag)
  {
    for (size_t i = 0; i < m.n_rows; i++)
    {
      for (size_t j = 0; j < i; j++)
      {
        std::cout << sep;
      }
      for (size_t j = i; j < m.n_cols; j++)
      {
        std::cout << m(i, j)
                  << ( j == m.n_cols - 1 ? '\n' : sep);
      }
    }
  }
  else
  {
    for (size_t i = 0; i < m.n_rows - 1; i++)
    {
      for (size_t j = 0; j < i; j++)
      {
        std::cout << sep;
      }
      for (size_t j = i + 1; j < m.n_cols; j++)
      {
        std::cout << m(i, j)
                  << ( j == m.n_cols - 1 ? '\n' : sep);
      }
    }
  }
}

void write_aracne(mat_t& m, std::vector<std::string>& g,
                  char sep)
{
  sep = (sep == '\0' ? '\t' : sep);
  for (size_t i = 0; i < m.n_rows; i++)
  {
    std::cout << g[i] << sep;
    for (size_t j = 0; j < m.n_cols; j++)
    {
      std::cout << g[j] << sep;
      std::cout << m(i, j)
                << (j == m.n_cols - 1 ? '\n' : sep);
    }
  }
}

int convert(int argc, char ** argv)
{
  LOG_INIT_CLOG();
  std::string infile;
  std::string gene_file;
  std::string fill;
  std::string in_format;
  std::string in_sep;
  std::string out_format;
  std::string out_sep;
  uint16_t prec;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " convert";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    tclap_seidr_format fmt;
    std::vector<std::string> formats{"edge-list", "sym-mat", "low-tri",
                                     "up-tri", "low-tri-diag", "up-tri-diag",
                                     "aracne"};
    TCLAP::ValuesConstraint<std::string> constraints_formats(formats);

    // Add arguments from the command line
    TCLAP::CmdLine cmd("Convert different text based formats\n"
                       "------------------------------------\n"
                       "This tool is meant to convert\n"
                       "between various text based network\n"
                       "formats. Currently supports parsing\n"
                       "and writing of:\n"
                       "  aracne\n"
                       "    The file format output by the\n"
                       "    ARACNE GRN software\n"
                       "  edge-list\n"
                       "    A three column file in the form\n"
                       "    <source> <target> <weight>\n"
                       "  low-tri\n"
                       "    A lower triangular matrix of\n"
                       "    scores without a diagonal\n"
                       "  low-tri-diag\n"
                       "    A lower triangular matrix of\n"
                       "    scores with a diagonal\n"
                       "  sym-mat\n"
                       "    A symmetric  matrix of scores\n"
                       "  up-tri\n"
                       "    An upper triangular matrix of\n"
                       "    scores without a diagonal\n"
                       "  up-tri-diag\n"
                       "    An upper triangular matrix of\n"
                       "    scores with a diagonal\n",
                       ' ', version);
    cmd.setOutput(&fmt);
    TCLAP::ValueArg<std::string>
    arg_infile("i", "infile", "Input file", true, "", "");

    TCLAP::ValueArg<std::string>
    arg_genefile("g", "genes", "Gene names", true, "", "");

    TCLAP::ValueArg<std::string>
    arg_fill("F", "fill", "Fill value for missing data (Number, NaN, Inf, -Inf)",
             false, "NaN", "NaN");

    TCLAP::ValueArg<std::string>
    arg_informat("f", "from", "Input file format", true, "", &constraints_formats);

    TCLAP::ValueArg<std::string>
    arg_outformat("t", "to", "Output file format", true, "", &constraints_formats);

    TCLAP::ValueArg<std::string>
    arg_outsep("S", "out-separator", "Output field separator",
               false, "\0", "tab");

    TCLAP::ValueArg<std::string>
    arg_insep("s", "in-separator", "Input field separator",
              false, "\0", "tab");

    TCLAP::ValueArg<uint16_t>
    arg_prec("p", "precision", "Number of significant digits to report", false,
             8, "8");

    cmd.add(arg_prec);
    cmd.add(arg_infile);
    cmd.add(arg_genefile);
    cmd.add(arg_fill);
    cmd.add(arg_informat);
    cmd.add(arg_outformat);
    cmd.add(arg_outsep);
    cmd.add(arg_insep);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    gene_file = arg_genefile.getValue();
    fill = arg_fill.getValue();
    in_format = arg_informat.getValue();
    out_format = arg_outformat.getValue();
    out_sep = arg_outsep.getValue();
    in_sep = arg_insep.getValue();
    prec = arg_prec.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log << except.what() << '\n';
    return EINVAL;
  }

  std::cout.precision(prec);
  std::cout.setf( std::ios::fixed, std:: ios::floatfield );

  std::vector<std::string> genes = read_genes(gene_file, '\n', '\t');
  seidr_score_t fi;
  if (fill == "Inf")
    fi = arma::datum::inf;
  else if (fill == "-Inf")
    fi = -arma::datum::inf;
  else if (fill == "NaN")
    fi = arma::datum::nan;
  else
  {
#ifdef SEIDR_SCORE_DOUBLE
    fi = std::stod(fill);
#else
    fi = std::stof(fill);
#endif
  }

  mat_t gm(genes.size(), genes.size());

  gm.fill(fi);

  std::map<std::string, size_t> gene_map;
  for (size_t i = 0; i < genes.size(); i++)
  {
    gene_map[ genes[i] ] = i;
  }

  std::shared_ptr<std::istream> in_stream;
  if (infile == "-")
    in_stream.reset(&(std::cin), [](...) {});
  else
    in_stream.reset(new std::ifstream(infile.c_str(), std::ios::in ) );

  // Generic reading/writing if no specialized function is available
  if (in_format == "edge-list")
    parse_el(gm, gene_map, *in_stream, in_sep[0]);

  else if (in_format == "sym-mat")
    parse_sm(gm, *in_stream, in_sep[0]);

  else if (in_format == "low-tri")
    parse_ltri(gm, *in_stream, false, in_sep[0]);

  else if (in_format == "up-tri")
    parse_utri(gm, *in_stream, false, in_sep[0]);

  else if (in_format == "low-tri-diag")
    parse_ltri(gm, *in_stream, true, in_sep[0]);

  else if (in_format == "up-tri-diag")
    parse_utri(gm, *in_stream, true, in_sep[0]);

  else if (in_format == "aracne")
    parse_aracne(gm, gene_map, *in_stream, in_sep[0]);


  // Output
  if (out_format == "edge-list")
    write_el(gm, genes, fi, out_sep[0]);

  else if (out_format == "sym-mat")
    write_sm(gm, out_sep[0]);

  else if (out_format == "low-tri")
    write_ltri(gm, false, out_sep[0]);

  else if (out_format == "up-tri")
    write_utri(gm, false, out_sep[0]);

  else if (out_format == "low-tri-diag")
    write_ltri(gm, true, out_sep[0]);

  else if (out_format == "up-tri-diag")
    write_utri(gm, true, out_sep[0]);

  else if (out_format == "aracne")
    write_aracne(gm, genes, out_sep[0]);

  return 0;
}
