/**
 * @file
 * @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
 * @version 0.01
 *
 * @section DESCRIPTION
 *
 * This is a collection of functions common to most other routines.
 */
// Seidr
#include <common.h>
#include <fs.h>
#include <mpims.h>
// External
#include <armadillo>
#include <cerrno>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <map>

namespace fs = boost::filesystem;

typedef std::pair<std::streampos, std::string> file_index_t;
typedef std::map<seidr_uword_t, file_index_t> result_file_map_t;

/**
 * Get a set difference of a single gene versus all genes
 *
 * This function will return the set difference between any
 * single gene and the collection of all genes. Useful for
 * all vs. all style tests.
 *
 * @param ind A vector (armadillo uword) for indices in 1:N
 *            for all genes.
 * @param s   A single index which should be excluded from
 *            the output.
 * @return arma::uvec A vector of indicies without the
 *                    excluded one.
 */
arma::uvec get_i(arma::uword ind, size_t s) {
  arma::uvec res(s - 1);
  arma::uword ri = 0;
  for (arma::uword f = 0; f < s; f++) {
    if (f < ind) res(ri++) = f;
    if (f > ind) res(ri++) = f;
  }
  return res;
}

/**
 * Read in genes from a genes file.
 *
 * This function will read all genes into a vector of strings.
 * The gene names can be delimited by a tab or by newlines.
 *
 * @param input The file path.
 * @param row_delim The dilimiter character for rows ('\n')
 * @param field_delim The delimiter character for columns ('\t')
 * @return std::vector<std::string> of gene names.
 */

std::vector<std::string> read_genes(std::string input,
                                    char row_delim, char field_delim) {
  std::vector<std::string> res;
  const char rd = row_delim;
  const char fd = field_delim;

  std::ifstream inf;
  inf.open(input.c_str());
  for (std::string row; std::getline(inf, row, rd); )
  {
    std::istringstream ss(row);

    for (std::string field; std::getline(ss, field, fd); )
    {
      res.push_back(field);
    }
  }

  inf.close();

  return res;
}

/**
 * Test if a file is GZipped based on magic numbers
 *
 * @param input An iostream (open in binary)
 * @return bool true if the file is GZipped, false if not
 */

bool is_gzip(std::string input)
{
  std::ifstream _input(input.c_str(), std::ios::in | std::ios::binary);
  bool gzip = false;
  char byte1 = _input.get();
  char byte2 = _input.get();
  if ((byte1 == '\x1F') && (byte2 == '\x8B')) gzip = true;
  _input.close();
  return gzip;
}


//Almost equal a la Bruce Dawson
bool almost_equal(seidr_score_t A, seidr_score_t B)
{
  // Calculate the difference.
  seidr_score_t diff = fabs(A - B);
  A = fabs(A);
  B = fabs(B);
  // Find the largest
  seidr_score_t largest = (B > A) ? B : A;

  if (diff <= largest * FLT_EPSILON)
    return true;
  return false;
}

seidr_score_t unity_stand(seidr_score_t xmin, seidr_score_t xmax, seidr_score_t xi)
{
  seidr_score_t x = 1 - ( (xi - xmin) / (xmax - xmin) );
  return x;
}

void scale(arma::mat& x) {
  x.each_col([](arma::vec & v) {
    double m = arma::mean(v);
    double sd = arma::stddev(v);
    v -= m;
    v /= sd;
  });
}

std::vector<std::string> tokenize_delim(std::string nodes, std::string delim)
{
  std::vector<std::string> nodelist;
  boost::char_separator<char> sep(delim.c_str());
  tokenizer tokens(nodes, sep);
  for (auto tok = tokens.begin(); tok != tokens.end(); tok++)
  {
    nodelist.push_back(*tok);
  }
  return nodelist;
}

void merge_files(std::string outfile, std::string outfilebase,
                 std::string tempdir, bool targeted, int id,
                 std::vector<std::string>& genes)
{
  result_file_map_t rmap;
  if (id == 0)
  {
    seidr_mpi_logger log;
    log << "Merging tmp files and cleaning up.\n";
    log.send(LOG_INFO);
    std::vector<fs::path> files;
    std::string tmpdir = outfilebase + "/" + tempdir;
    fs::path p_tmp(tmpdir);
    for (auto it = fs::directory_iterator(p_tmp);
         it != fs::directory_iterator(); it++)
    {
      if ( fs::is_regular_file( it->path() ) )
        files.push_back( (*it).path() );
    }
    std::ofstream ofs(outfile);

    for (fs::path& p : files)
    {
      std::ifstream ifs(p.string().c_str());
      std::string l;
      while (std::getline(ifs, l))
      {
        seidr_uword_t gene_index = std::stoul(l);
        std::streampos g = ifs.tellg();
        rmap[gene_index] = file_index_t(g, p.string());
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      ifs.close();
    }

    auto it = rmap.begin();
    auto gene_index = it->second.first;
    auto file_path = it->second.second;
    std::ifstream ifs(file_path);
    for (; it != rmap.end(); it++)
    {
      ifs.seekg(gene_index);
      std::string l;
      std::getline(ifs, l);

      if (targeted)
      {
        std::stringstream ss(l);
        std::string token;
        seidr_uword_t j = 0;
        seidr_uword_t i = it->first;
        while (ss >> token)
        {
          if (i != j)
            ofs << genes[i] << '\t' << genes[j] << '\t' << token << '\n';
          j++;
        }
      }
      else
      {
        ofs << l << '\n';
      }

      auto nx = std::next(it);
      if (nx != rmap.end())
      {
        if (nx->second.second != file_path)
        {
          ifs.close();
          ifs.open(nx->second.second);
          gene_index = nx->second.first;
        }
        else
        {
          gene_index = nx->second.first;
        }
      }
    }
    remove(tmpdir, true);
  }
}

bool any_const_expr(arma::mat& inp)
{
  arma::mat v = arma::var(inp);
  for (arma::uword i = 0; i < v.n_elem; i++)
  {
    if (almost_equal(v(0, i), 0))
    {
      return true;
    }
  }
  return false;
}

void verify_matrix(arma::mat& inp)
{
  if (any_const_expr(inp))
  {
    throw std::runtime_error("Constant values detected in at least one column"
                               ". Please filter your input to contain only "
                               "columns with non-zero variance.");
  }
}