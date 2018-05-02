//Seidr
#include <mi_fun.h>
#include <mpims.h>
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

using boost::lexical_cast;
namespace fs = boost::filesystem;

bool sortDesc(const aranode& lhs, const aranode& rhs)
{
  return lhs.v > rhs.v;
}

class seidr_mpi_mi : public seidr_mpi
{
public:
  void entrypoint();
  void finalize();
  seidr_mpi_mi() : seidr_mpi() {}
  seidr_mpi_mi(unsigned long bs) : seidr_mpi(bs) {}
  void set_num_bins(size_t n) {_num_bins = n;}
  void set_spline_order(size_t s) {_spline_order = s;}
  void set_mode(char x) {_mode = x;}
  void set_mi_file(std::string x) {_mi_file = x;}
private:
  size_t _num_bins;
  size_t _spline_order;
  char _mode;
  std::string _tmpdir;
  std::string _mi_file;
};

void seidr_mpi_mi::entrypoint()
{
  seidr_mpi_logger log;
  std::string mode;
  switch (_mode)
  {
  case 0:
    mode = std::string("CLR");
    break;
  case 1:
    mode = std::string("ARACNE");
    break;
  case 2:
    mode = std::string("RAW");
    break;
  default:
    throw std::runtime_error("Post processing mode not known");
  }
  _tmpdir = _outfilebase + "/.seidr_tmp_mi_" + mode;
  std::string tmpfile = _tmpdir + "/MPIthread_" +
                        std::to_string(_id) + "_" +
                        std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);

  std::vector<arma::uword> uvec;
  for (auto i : _indices)
    uvec.push_back(i);
  mi_sub_matrix(_data, _num_bins, _spline_order, uvec, tmpfile);
  announce_ready();
}

void seidr_mpi_mi::finalize()
{
  if (_id == 0)
  {
    seidr_mpi_logger log;
    log << "Merging tmp files\n";
    log.log(LOG_INFO);
    std::vector<fs::path> files;
    fs::path p_tmp(_tmpdir);
    for (auto it = fs::directory_iterator(p_tmp);
         it != fs::directory_iterator(); it++)
    {
      if ( fs::is_regular_file( it->path() ) )
        files.push_back( (*it).path() );
    }
    std::ofstream ofs(_outfile);

    arma::mat mi_mat(_data.n_cols, _data.n_cols, arma::fill::zeros);

    for (fs::path& p : files)
    {
      std::ifstream ifs(p.string().c_str());
      std::string l;
      while (std::getline(ifs, l))
      {
        size_t row = std::stoul(l);
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

    if (_mode == 2 || _mi_file != "")
    {
      log << "Outputting raw MI\n";
      log.send(LOG_INFO);
      bool only_mi = true;
      std::ofstream mofs;
      if (_mi_file != "")
      {
        only_mi = false;
        mofs.open(_mi_file.c_str());
      }
      for (arma::uword i = 1; i < mi_mat.n_cols; i++)
      {
        for (arma::uword j = 0; j < i; j++)
        {

          (only_mi ? ofs : mofs) << mi_mat(i, j)
                                 << (j == i - 1 ? '\n' : '\t');
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
    log << "Finished\n";
    log.log(LOG_INFO);
    fs::remove_all(_tmpdir);
  }
  // get remaining logs
  double now = MPI_Wtime();
  while (MPI_Wtime() - now < 10)
  {
    check_logs();
  }
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

arma::uvec bin_count(arma::mat& gm, size_t multiplier)
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


double basis_function(size_t i, size_t p, double t, arma::vec& knot_vector, size_t num_bins)
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

void find_weights(arma::mat& gm, arma::vec& knots, arma::mat& wm,
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

double entropy1d(arma::mat& gm, arma::vec& knots, arma::mat& wm,
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
          wy(cur_bin_y * x.n_elem + cur_sample) /
          n;
      }
    }
  }

}

double entropy2d(arma::mat& gm, arma::vec& knots,
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

void mi_sub_matrix(arma::mat& gm, size_t num_bins, size_t spline_order,
                   std::vector<arma::uword> targets, std::string tempfile)
{
  seidr_mpi_logger log;
  arma::vec knots = knot_vector(spline_order, num_bins);
  arma::mat weights(gm.n_rows * num_bins, gm.n_cols, arma::fill::zeros);
  arma::vec entropies(gm.n_cols, arma::fill::zeros);

  std::sort(targets.begin(), targets.end());

  std::ofstream ofs(tempfile.c_str(), std::ios::out);

  for (size_t i = 0; i < gm.n_cols; i++)
  {
    find_weights(gm, knots, weights, spline_order, num_bins, i);
    entropies(i) = entropy1d(gm, knots, weights, spline_order, num_bins, i);
  }

  for (size_t row : targets)
  {
    if (row > 0)
    {
      ofs << row << '\n';
      for (size_t col = 0; col < row; col++)
      {
        double e2d = entropy2d(gm, knots, weights, spline_order, num_bins,
                               col, row);
        ofs <<  entropies(col) + entropies(row) - e2d;
        ofs << (col == row - 1 ? '\n' : '\t');
      }
      log << "Finished gene " << row << '\n';
      log.send(LOG_INFO);
    }
  }
}


void mi_full(arma::mat& gm, size_t spline_order, size_t num_bins, size_t bs,
             std::string outfile, char mode, std::string mi_file)
{
  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<unsigned long> targets(gm.n_cols);
  for (size_t i = 0; i < gm.n_cols; i++)
    targets[i] = i;

  std::random_shuffle(targets.begin(), targets.end());

  seidr_mpi_mi mpi(bs);

  mpi.set_data(gm);
  mpi.set_indices(targets);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_spline_order(spline_order);
  mpi.set_num_bins(num_bins);
  mpi.set_mode(mode);
  mpi.set_mi_file(mi_file);

  mpi.exec();
  //mpi.finalize();
}
