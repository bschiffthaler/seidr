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

//Seidr
#include <mi_fun.h>
#include <mpiomp.h>
//External
#include <armadillo>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>

using boost::lexical_cast;
namespace fs = boost::filesystem;

bool sortDesc(const aranode& lhs, const aranode& rhs)
{
  return lhs.v > rhs.v;
}

class seidr_mpi_mi : public seidr_mpi_omp
{
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_num_bins(size_t n) {_num_bins = n;}
  void set_spline_order(size_t s) {_spline_order = s;}
  void set_mode(char x) {_mode = x;}
  void set_mi_file(std::string x) {_mi_file = x;}
  void set_targets(std::vector<std::string> x) {_targets = x;}
  void set_genes(std::vector<std::string> x) {_genes = x;}
  void use_existing_mi_mat() {_use_existing_mi_mat = true;}
private:
  size_t _num_bins = 0;
  size_t _spline_order = 0;
  char _mode = 0;
  bool _use_existing_mi_mat = false;
  std::string _mi_file = "";
  std::vector<std::string> _targets;
  std::vector<std::string> _genes;
};

void seidr_mpi_mi::entrypoint()
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  if (! _use_existing_mi_mat)
  {
    while (! _my_indices.empty())
    {
      std::vector<arma::uword> uvec;
      for (auto i : _my_indices)
      {
        uvec.push_back(i);
      }
      mi_sub_matrix(_data, _num_bins, _spline_order, uvec, _tempdir);
      get_more_work();
    }
    #pragma omp critical
    {
      log << "No more work. Waiting for other tasks to finish...\n";
      log.send(LOG_INFO);
    }
  }
}

void seidr_mpi_mi::finalize()
{
  arma::mat mi_mat(_data.n_cols, _data.n_cols, arma::fill::zeros);
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  std::unordered_map<std::string, arma::uword> gene_map;
  arma::uword ctr = 0;
  for (auto g : _genes)
  {
    gene_map[g] = ctr++;
  }

  if (_use_existing_mi_mat)
  {
    log << "Using MI from " << _mi_file << '\n';
    log.log(LOG_INFO);
    std::ifstream mifs(_mi_file.c_str());
    std::string l;
    arma::uword i = 1;
    while (std::getline(mifs, l))
    {
      arma::uword j = 0;
      std::string f;
      std::stringstream ss(l);
      while (std::getline(ss, f, '\t'))
      {
        double val = std::stod(f);
        mi_mat(i, j) = val;
        mi_mat(j, i) = val;
        j++;
      }
      i++;
    }
  }
  else
  {
    log << "Merging tmp files from " << _tempdir << '\n';
    log.log(LOG_INFO);

    std::vector<fs::path> files;
    fs::path p_tmp(_tempdir);

    for (auto it = fs::directory_iterator(p_tmp);
         it != fs::directory_iterator(); it++)
    {
      std::string pstring = it->path().string();
      if (pstring.find(".") != std::string::npos)
      {
        log << "Ignoring unexpected file: "
            << pstring << '\n';
        log.log(LOG_WARN);
        while (check_logs(LOG_NAME"@" + mpi_get_host()));
      }
      else if ( fs::is_regular_file( it->path() ) )
      {
        files.push_back( (*it).path() );
      }
    }

    for (fs::path& p : files)
    {
      std::ifstream ifs(p.string().c_str());
      std::string l;
      while (std::getline(ifs, l))
      {
        size_t row;
        try
        {
          row = std::stoul(l);
        }
        catch (std::exception& e)
        {
          throw e;
        }
        std::getline(ifs, l);
        std::stringstream ss(l);
        std::string field;
        size_t col = 0;
        while (std::getline(ss, field, '\t'))
        {
          mi_mat(row, col) = std::stod(field);
          mi_mat(col, row) = mi_mat(row, col);
          col++;
        }
      }
      log << "Merged " << p.string() << '\n';
      log.log(LOG_INFO);
    }

    log << "Merged all files.\n";
    log.send(LOG_INFO);
  }

  std::ofstream ofs(_outfile);

  if (_mi_file != "" && (! _use_existing_mi_mat))
  {
    std::ofstream mofs;
    mofs.open(_mi_file.c_str());
    log << "Outputting raw MI\n";
    log.send(LOG_INFO);
    for (arma::uword i = 1; i < mi_mat.n_cols; i++)
    {
      for (arma::uword j = 0; j < i; j++)
      {
        mofs << mi_mat(i, j) << (j == i - 1 ? '\n' : '\t');
      }
    }
    mofs.close();
  }

  if (_mode == 2)
  {
    if (_targets.size() != 0)
    {
      for (auto g : _targets)
      {
        arma::uword i;
        try
        {
          i = gene_map.at(g);
        }
        catch (std::exception& e)
        {
          log << e.what() << '\n';
          log << "Target gene " << g << " is not in expression table\n";
          log.send(LOG_ERR);
        }
        for (arma::uword j = 0; j < mi_mat.n_cols; j++)
        {
          if (i == j) continue;
          ofs << _genes[i] << '\t'
              << _genes[j] << '\t'
              << mi_mat(i, j) << '\n';
        }
      }
    }
    else
    {
      for (arma::uword i = 1; i < mi_mat.n_cols; i++)
      {
        for (arma::uword j = 0; j < i; j++)
        {
          ofs << mi_mat(i, j) << (j == i - 1 ? '\n' : '\t');
        }
      }
    }
  }
  if (_mode == 0)
  {
    log << "Computing CLR\n";
    log.log(LOG_INFO);

    arma::vec m(mi_mat.n_cols, arma::fill::zeros);
    arma::vec s(mi_mat.n_cols, arma::fill::zeros);

    for (size_t i = 0; i < mi_mat.n_cols; i++)
    {
      arma::vec x = mi_mat.col(i);
      x = arma::abs(x);
      m(i) = arma::mean(x);
      s(i) = arma::stddev(x);
    }

    if (_targets.size() != 0)
    {
      for (auto g : _targets)
      {
        arma::uword i;
        try
        {
          i = gene_map.at(g);
        }
        catch (std::exception& e)
        {
          log << e.what() << '\n';
          log << "Target gene " << g << " is not in expression table\n";
          log.send(LOG_ERR);
        }
        for (arma::uword j = 0; j < mi_mat.n_cols; j++)
        {
          if (i == j) continue;
          double mi = mi_mat(i, j);
          mi = mi < 0 ? mi * -1 : mi;
          double a = (mi - m(i)) / s(i);
          double b = (mi - m(j)) / s(j);
          a = a < 0 ? 0 : a;
          b = b < 0 ? 0 : b;
          double c = sqrt(a * a + b * b);
          c = c < 0 ? 0 : c;
          ofs
              << _genes[i] << '\t'
              << _genes[j] << '\t'
              << c << '\n';
        }
      }
    }
    else
    {
      for (size_t i = 1; i < mi_mat.n_cols; i++)
      {
        for (size_t j = 0; j < i; j++)
        {
          double mi = mi_mat(i, j);
          mi = mi < 0 ? mi * -1 : mi;
          double a = (mi - m(i)) / s(i);
          double b = (mi - m(j)) / s(j);
          a = a < 0 ? 0 : a;
          b = b < 0 ? 0 : b;
          double c = sqrt(a * a + b * b);
          c = c < 0 ? 0 : c;
          ofs << c;
          ofs << (j == i - 1 ? '\n' : '\t');
        }
      }
    }
  }
  else if (_mode == 1)
  {
    log << "Applying DPI to data\n";
    log.send(LOG_INFO);
    arma::imat torm(mi_mat.n_cols, mi_mat.n_cols, arma::fill::zeros);
    for (arma::uword i = 0; i < mi_mat.n_rows; i++)
    {
      std::vector<aranode> dc;
      for (arma::uword j = 0; j < mi_mat.n_cols; j++)
      {
        if (i != j)
          dc.push_back(aranode(j, mi_mat(i, j)));
      }
      std::sort(dc.begin(), dc.end(), sortDesc);
      for (arma::uword j = 0; j < dc.size(); j++)
      {
        double ab = dc[j].v;
        arma::uword bi = dc[j].i;
        for (arma::uword k = 0; k < j; k++)
        {
          double ac = dc[k].v;
          if (ab >= ac)
            break;
          arma::uword ci = dc[k].i;
          double bc = mi_mat(bi, ci);
          if (ab < bc)
          {
            torm(i, bi) = 1;
            torm(bi, i) = 1;
            break;
          }
        }
      }
      if (i > 0 && i % 100 == 0)
      {
        log << "Processed " << i << '/' << mi_mat.n_rows << '\n';
        log.send(LOG_INFO);
      }
    }
    if (_targets.size() != 0)
    {
      for (auto g : _targets)
      {
        arma::uword i;
        try
        {
          i = gene_map.at(g);
        }
        catch (std::exception& e)
        {
          log << e.what() << '\n';
          log << "Target gene " << g << " is not in expression table\n";
          log.send(LOG_ERR);
        }
        for (arma::uword j = 0; j < mi_mat.n_cols; j++)
        {
          if (i == j) continue;
          double r = (torm(i, j) ? 0 : mi_mat(i, j));
          ofs
              << _genes[i] << '\t'
              << _genes[j] << '\t'
              << r << '\n';
        }
      }
    }
    else
    {
      for (arma::uword i = 1; i < mi_mat.n_cols; i++)
      {
        for (arma::uword j = 0; j < i; j++)
        {
          double r = (torm(i, j) ? 0 : mi_mat(i, j));
          ofs << r
              << (j == i - 1 ? '\n' : '\t');
        }
      }
    }
  }
  log << "Finished\n";
  log.log(LOG_INFO);
  fs::remove_all(_tempdir);
}


arma::vec knot_vector(size_t spline_order, size_t num_bins)
{
  arma::vec ret(spline_order + num_bins, arma::fill::zeros);
  size_t int_points = num_bins - spline_order;

  for (size_t i = spline_order; i < spline_order + int_points; i++)
  {
    double x = (i - spline_order + 1);
    double y = (int_points + 1);
    ret[i] = x / y;
  }

  for (size_t i = spline_order + int_points; i < 2 * spline_order + int_points; i++)
  {
    ret[i] = 1;
  }

  return ret;
}

/* Linear interpolation percentile method as implemented
   in e.g. MATLAB. Requires data to be SORTED
*/
double percentile(arma::vec& data, size_t percentile)
{
  std::vector<size_t> p_rank{};
  double N = data.n_elem;

  for (size_t i = 0; i < data.n_elem; i++)
  {
    double di = (i + 1);
    p_rank.push_back( round((100 / N) * (di - 0.5) ));
  }

  // percentile is in ordered list
  if (percentile <= p_rank[0]) return data[0];
  if (percentile >= p_rank[p_rank.size() - 1])
    return  data[p_rank.size() - 1];

  auto it = std::lower_bound(p_rank.begin(), p_rank.end(), percentile);

  if ( (*it) == percentile )
    return data[std::distance(p_rank.begin(), it)];

  // percentile is not in ordered list
  size_t k2   = std::distance(p_rank.begin(), it);
  size_t k    = std::distance(p_rank.begin(), --it);
  double pk   = p_rank[k];
  double vk   = data[k];
  double vk2  = data[k2];
  return vk + N * ((percentile - pk) / 100) * (vk2 - vk);
}

/* Compute the interquartile range of a SORTED
   arma::vec
*/
double iqr(arma::vec& data)
{
  return percentile(data, 75) - percentile(data, 25);
}

double bin_width(arma::vec& data)
{
  double s = arma::stddev(data);
  double i = iqr(data) / 1.349;
  double sigma = s < i ? s : i;
  double ns = data.n_elem;
  return pow(ns, (-1.0 / 3.0)) * sigma * 3.49;
}

arma::uvec bin_count(const arma::mat& gm, size_t multiplier)
{
  arma::uvec ret(gm.n_cols, arma::fill::zeros);

  arma::mat s_data = gm;
  double m = multiplier;

  for (size_t i = 0; i < gm.n_cols; i++)
  {
    arma::vec v = gm.col(i);
    v = arma::sort(v);
    double x = v[v.size() - 1] - v[0];
    double y = bin_width(v) * m;
    ret(i) = ceil(x / y);
  }

  ret = arma::sort(ret);

  return ret;
}

arma::vec to_z(arma::vec& x) {
  double mi = arma::min(x);
  double ma = arma::max(x);
  arma::vec z = x;
  for (arma::vec::iterator i = z.begin(); i != z.end(); i++)
    (*i) = ( ( (*i) - mi ) / ( ma - mi ) );
  return z;
}


double basis_function(size_t i, size_t p, double t, arma::vec& knot_vector,
                      size_t num_bins)
{
  double d1, n1, d2, n2, e1, e2;
  if (p == 1)
  {
    if ((t >= knot_vector[i] && t < knot_vector[i + 1] &&
         knot_vector[i] < knot_vector[i + 1]) ||
        (fabs(t - knot_vector[i + 1]) < 1e-10 && (i + 1 == num_bins)))
    {
      return (1);
    }
    return (0);
  }

  d1 = knot_vector[i + p - 1] - knot_vector[i];
  n1 = t - knot_vector[i];
  d2 = knot_vector[i + p] - knot_vector[i + 1];
  n2 = knot_vector[i + p] - t;

  if (d1 < 1e-10 && d2 < 1e-10)
  {
    return (0);
  }
  else if (d1 < 1e-10)
  {
    e1 = 0;
    e2 = n2 / d2 * basis_function(i + 1, p - 1, t, knot_vector, num_bins);
  }
  else if (d2 < 1e-10)
  {
    e2 = 0;
    e1 = n1 / d1 * basis_function(i, p - 1, t, knot_vector, num_bins);
  }
  else
  {
    e1 = n1 / d1 * basis_function(i, p - 1, t, knot_vector, num_bins);
    e2 = n2 / d2 * basis_function(i + 1, p - 1, t, knot_vector, num_bins);
  }

  /* sometimes, this value is < 0 (only just; rounding error); truncate */
  if (e1 + e2 < 0)
  {
    return (0);
  }
  return (e1 + e2);

}

void find_weights(const arma::mat& gm, arma::vec& knots, arma::mat& wm,
                  size_t spline_order, size_t num_bins, size_t i)
{
  // standardize data
  arma::vec z = gm.col(i);
  z = to_z(z);
  // insert weights for each bin into the weight matrix
  for (size_t cur_sample = 0; cur_sample < z.n_elem; cur_sample++)
  {
    for (size_t cur_bin = 0; cur_bin < num_bins; cur_bin++)
    {
      wm(cur_bin * z.n_elem + cur_sample, i) =
        basis_function(cur_bin, spline_order, z[cur_sample], knots, num_bins);
    }
  }
}


arma::vec hist1d(arma::vec& x, arma::vec& knots, arma::vec& weights,
                 size_t spline_order, size_t num_bins)
{
  arma::vec hist(num_bins, arma::fill::zeros);
  for (size_t cur_bin = 0; cur_bin < num_bins; cur_bin++)
  {
    for (size_t cur_sample = 0; cur_sample < x.n_elem; cur_sample++)
    {
      double n = x.n_elem;
      hist[cur_bin] += weights[cur_bin * x.n_elem + cur_sample] / n;
    }
  }
  return hist;
}

double log2d(double x)
{
  return log(x) / log(2);
}

double entropy1d(const arma::mat& gm, arma::vec& knots, arma::mat& wm,
                 size_t spline_order, size_t num_bins, size_t i)
{
  double H = 0;
  arma::vec weights = wm.col(i);
  arma::vec x = gm.col(i);
  arma::vec hist = hist1d(x, knots, weights, spline_order, num_bins);
  for (size_t cur_bin = 0; cur_bin < num_bins; cur_bin++)
  {
    if (hist[cur_bin] > 0)
    {
      H -= hist[cur_bin] * log2d(hist[cur_bin]);
    }
  }
  return H;
}

void hist2d(arma::vec& x, arma::vec& y, arma::vec& knots,
            arma::vec& wx, arma::vec& wy, arma::mat& hist,
            size_t spline_order, size_t num_bins)
{

  for (size_t cur_bin_x = 0; cur_bin_x < num_bins; cur_bin_x++)
  {
    for (size_t cur_bin_y = 0; cur_bin_y < num_bins; cur_bin_y++)
    {
      for (size_t cur_sample = 0; cur_sample < x.n_elem; cur_sample++)
      {
        double n = x.n_elem;
        hist(cur_bin_x, cur_bin_y) +=
          wx(cur_bin_x * x.n_elem + cur_sample) *
          wy(cur_bin_y * x.n_elem + cur_sample) / n;
      }
    }
  }

}

double entropy2d(const arma::mat& gm, arma::vec& knots,
                 arma::mat& wm, size_t spline_order,
                 size_t num_bins, size_t xi, size_t yi)
{

  arma::vec x = gm.col(xi);
  arma::vec y = gm.col(yi);
  arma::vec wx = wm.col(xi);
  arma::vec wy = wm.col(yi);
  arma::mat hist(num_bins, num_bins, arma::fill::zeros);
  double H = 0;
  double incr;
  hist2d(x, y, knots, wx, wy, hist, spline_order, num_bins);

  for (size_t cur_bin_x = 0; cur_bin_x < num_bins; cur_bin_x++)
  {
    for (size_t cur_bin_y = 0; cur_bin_y < num_bins; cur_bin_y++)
    {
      incr = hist(cur_bin_x, cur_bin_y);
      if (incr > 0) {
        H -= incr * log2d(incr);
      }
    }
  }

  return H;
}

void mi_sub_matrix(const arma::mat& gm, size_t num_bins, size_t spline_order,
                   std::vector<arma::uword>& targets,
                   const std::string& tmpdir)
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  arma::vec knots = knot_vector(spline_order, num_bins);
  arma::mat weights(gm.n_rows * num_bins, gm.n_cols, arma::fill::zeros);
  arma::vec entropies(gm.n_cols, arma::fill::zeros);

  std::sort(targets.begin(), targets.end());

  #pragma omp parallel for
  for (size_t i = 0; i < gm.n_cols; i++)
  {
    find_weights(gm, knots, weights, spline_order, num_bins, i);
    entropies(i) = entropy1d(gm, knots, weights, spline_order, num_bins, i);
  }

  std::string tmpfile = tempfile(tmpdir);
  std::ofstream ofs(tmpfile.c_str(), std::ios::out);

  #pragma omp parallel for
  for (size_t i = 0; i < targets.size(); i++)
  {
    size_t row = targets[i];
    if (row > 0)
    {
      // Parallelize inner loop which works a bit better with pragma critical
      // and won't make much of a difference in practice as row == col
      std::vector<double> tmp_results;
      for (size_t col = 0; col < row; col++)
      {
        double e2d = entropy2d(gm, knots, weights, spline_order, num_bins,
                               col, row);
        tmp_results.push_back(entropies(col) + entropies(row) - e2d);
      }
      #pragma omp critical
      {
        ofs << row << '\n';
        for (uint64_t col = 0; col < row; col++)
        {
          ofs << tmp_results[col] << (col == row - 1 ? '\n' : '\t');
        }
        log << "Finished gene " << row << '\n';
        log.send(LOG_INFO);
      }
    }
  }
  ofs.close();
}


void mi_full(const arma::mat & gm,
             const std::vector<std::string>& genes,
             std::vector<std::string>& targets,
             const seidr_mi_param_t& param)
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  std::string mode_str;
  switch (param.m)
  {
  case 0:
    mode_str = std::string("CLR");
    break;
  case 1:
    mode_str = std::string("ARACNE");
    break;
  case 2:
    mode_str = std::string("RAW");
    break;
  default:
    throw std::runtime_error("Post processing mode not known");
  }

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> indices(gm.n_cols);
  for (uint64_t i = 0; i < gm.n_cols; i++)
  {
    indices[i] = i;
  }

  std::random_shuffle(indices.begin(), indices.end());

  seidr_mpi_mi mpi(param.bs, gm, indices, genes, param.tempdir,
                   param.outfile);
  mpi.set_spline_order(param.spline_order);
  mpi.set_num_bins(param.num_bins);
  mpi.set_mode(param.m);
  mpi.set_mi_file(param.mi_file);
  mpi.set_genes(genes);
  if (param.use_existing)
  {
    mpi.use_existing_mi_mat();
  }
  mpi.set_targets(targets);
  mpi.entrypoint();

  MPI_Barrier(MPI_COMM_WORLD); // NOLINT
  mpi.remove_queue_file();
  #pragma omp critical
  {
    if (mpi.rank() == 0)
    {
      while (mpi.check_logs(LOG_NAME"@" + mpi_get_host())); // NOLINT
      log << "Finalizing...\n";
      log.send(LOG_INFO);
      mpi.finalize();
      while (mpi.check_logs(LOG_NAME"@" + mpi_get_host()));
    }
  }

  MPI_Finalize();
}
