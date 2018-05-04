// Seidr
#include <narromi_fun.h>
#include <common.h>
#include <IP_LPT.h>
#include <stats_fun.h>
#include <mpims.h>
// External
#include <mpi.h>
#include <armadillo>
#include <vector>
#include <fstream>
#include <sys/types.h>
#include <sys/wait.h>
#include <ctime>
#include <boost/filesystem.hpp>
#include <map>

typedef std::pair<std::streampos, std::string> file_index_t;
typedef std::map<seidr_uword_t, file_index_t> result_file_map_t;

namespace fs = boost::filesystem;


namespace fs = boost::filesystem;

NaResult::NaResult(arma::vec sig, arma::vec nv,
                   arma::uvec re, arma::uword i) {
  s = sig; n = nv; r = re; I = i;
}
arma::vec NaResult::sig() { return s; }
arma::vec NaResult::nv() { return n; }
arma::uvec NaResult::re() { return r; }
arma::uword NaResult::i() { return I; }


class seidr_mpi_narromi : public seidr_mpi {
public:
  void entrypoint();
  void finalize();
  seidr_mpi_narromi() : seidr_mpi() {}
  seidr_mpi_narromi(unsigned long bs) : seidr_mpi(bs) {}
  void set_algorithm(std::string a) {_algorithm = a;}
  void set_alpha(double a) {_alpha = a;}
  void set_t(double t) {_t = t;}
  void set_targeted(bool x) {_targeted = x;}
private:
  std::string _algorithm;
  double _alpha;
  double _t;
  bool _targeted;
};

void seidr_mpi_narromi::entrypoint()
{
  seidr_mpi_logger log;
  std::string tmpfile = _outfilebase + "/.seidr_tmp_narromi/MPIthread_" +
                        std::to_string(_id) + "_" + std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);

  std::vector<arma::uword> uvec;
  for (auto i : _indices)
    uvec.push_back(i);
  narromi_thread(_data, _algorithm, _alpha, _t, uvec, _genes, tmpfile, _id);
  announce_ready();
}

void seidr_mpi_narromi::finalize()
{
  result_file_map_t rmap;
  if (_id == 0)
  {
    seidr_mpi_logger log;
    log << "Merging tmp files and cleaning up.\n";
    log.send(LOG_INFO);
    std::vector<fs::path> files;
    std::string tmpdir = _outfilebase + "/.seidr_tmp_narromi";
    fs::path p_tmp(tmpdir);
    for (auto it = fs::directory_iterator(p_tmp);
         it != fs::directory_iterator(); it++)
    {
      if ( fs::is_regular_file( it->path() ) )
        files.push_back( (*it).path() );
    }
    std::ofstream ofs(_outfile);

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

      if (_targeted)
      {
        std::stringstream ss(l);
        std::string token;
        seidr_uword_t j = 0;
        seidr_uword_t i = it->first;
        while(ss >> token)
        {
          if(i != j)
            ofs << _genes[i] << '\t' << _genes[j] << '\t' << token << '\n';
          j++;
        }
      }
      else
      {
        ofs << l << '\n';
      }

      auto nx = std::next(it);
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
    fs::remove_all(p_tmp);
    double now = MPI_Wtime();
    // Push all pending logs
    log << "Waiting up to 10 seconds for queued logs...\n";
    log.send(LOG_INFO);
    while (MPI_Wtime() - now < 10)
    {
      check_logs();
    }
  }
}

/* Linear programming optimisation loop.
 */
std::pair<arma::vec, arma::vec>
reoptim(arma::mat X, arma::rowvec Y,
        std::string al, double alpha) {

  arma::mat J = lp_ipt(X, Y, al);
  arma::vec Js(J.n_cols);
  arma::vec Jv(J.n_cols);

  Js.zeros(); Jv.zeros();

  arma::uvec index   = arma::find(arma::abs(J) >= alpha);
  arma::uvec index_c = arma::find(arma::abs(J) < alpha);

  Js(index) = J(index);
  Jv(index_c) = J(index_c);

  while (index_c.n_elem > 0) {
    arma::mat J1 = lp_ipt(X.rows(index), Y, al);
    arma::uvec index1   = arma::find(arma::abs(J1) >= alpha);
    arma::uvec index1_c = arma::find(arma::abs(J1) < alpha);
    index_c     = index(index1_c);
    index       = index(index1);
    Js(index)   = J1(index1);
    Js(index_c).zeros();
    Jv(index_c) = J1(index1_c);
  }

  Jv = Js + Jv;

  return std::pair<arma::vec, arma::vec> (Js, Jv);
}

/* Significance calculation function */

arma::vec get_sig(arma::vec z) {
  z = (z - arma::min(z)) / (arma::max(z) - arma::min(z));
  z = 0.5 * log2((1 + z) / (1 - z));
  // Set max values to second highest values
  arma::uvec m1 = arma::find(z == arma::max(z));
  arma::uvec m2 = arma::find(z != arma::max(z));
  z(m1).fill(arma::max(z(m2)));
  double mu = mean(z);
  double sigma = var(z);
  arma::vec sig = 1 - normcdf(z, mu, sigma);
  return sig;
}

// Run narromi algorithm for one gene against a set
// of genes
NaResult narromi(arma::mat GM, std::string al, double alpha, double t,
                 arma::uword i) {

  arma::uvec response = get_i(i, GM.n_cols);
  arma::mat X = GM.cols(response);
  arma::vec Y = GM.col(i);
  arma::vec net(X.n_cols, arma::fill::zeros);
  arma::vec fmi       = fullMI(X, Y);
  arma::vec net_value = fmi;
  arma::uvec index    = find(arma::abs(fmi) >= alpha);
  arma::mat X1        = X.cols(index);

  std::pair<arma::vec, arma::vec> J = reoptim(X1.t(), Y.t(), al, alpha);

  for (size_t i = 0; i < index.n_elem; i++) {
    net(index)        = J.first;
    net_value(index)  = J.second;
  }

  arma::vec net_value1 = net_value;
  net_value = arma::sign(net_value) % (arma::abs(net_value) * t + fmi * (1 - t));
  arma::vec z = arma::abs(net_value1);
  arma::vec sig1 = get_sig(z);

  z = (arma::abs(fmi));
  arma::vec sig2 = get_sig(z);

  arma::vec sig = arma::sqrt(arma::pow(sig1, 2) + arma::pow(sig2, 2));

  return NaResult(sig, net_value, response, i);

}

/* Run narromi algorithm on a full expression set in an all vs.
all style comparison.
 */
void full_narromi(arma::mat GM,
                  std::string al, double alpha, double t,
                  std::vector<std::string> genes, size_t bs,
                  std::string outfile) {

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<unsigned long> uvec;
  for (unsigned long i = 0; i < GM.n_cols; i++)
    uvec.push_back(i);

  seidr_mpi_narromi mpi(bs);

  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(uvec);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_alpha(alpha);
  mpi.set_t(t);
  mpi.set_algorithm(al);

  mpi.exec();
}

/* Run narromi algorithm on a full expression set on a set of
genes of interest.
*/

void partial_narromi(arma::mat GM,
                     std::string al, double alpha, double t,
                     std::vector<std::string> genes, size_t bs,
                     std::vector<std::string> targets,
                     std::string outfile) {
  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  seidr_mpi_logger log;

  std::vector<size_t> positions;
  for (size_t i = 0; i < targets.size(); i++) {
    size_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size()) {
      log << "\nWarning: Gene " << targets[i] << " was not found in the expression set " <<
          "and will therefore not be considered. Please check that your expression set and " <<
          "its column (gene) names contain an entry for " << targets[i] << ".\n";
      log.log(LOG_WARN);
    } else {
      positions.push_back(pos);
    }
  }

  seidr_mpi_narromi mpi(bs);

  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(positions);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_alpha(alpha);
  mpi.set_t(t);
  mpi.set_algorithm(al);
  mpi.set_targeted(true);

  mpi.exec();
}

void narromi_thread(arma::mat gene_matrix,
                    std::string al, double alpha, double t,
                    std::vector<arma::uword> ind,
                    std::vector<std::string> genes,
                    std::string tempfile, int id) {
  std::clock_t start;
  double time;
  size_t tot = gene_matrix.n_cols - 1;
  seidr_mpi_logger log;
  std::ofstream ofs(tempfile.c_str());

  if (! ofs)
    throw std::runtime_error("Could not open temp file: " + tempfile);

  for (size_t x = 0; x < ind.size(); x++) {
    start = std::clock();
    arma::uword i = ind[x];
    log << "Started gene " << genes[i] << '\n';
    log.send(LOG_INFO);
    NaResult nr = narromi(gene_matrix, al, alpha, t, i);
    arma::uvec re = nr.re();
    arma::vec si = nr.sig();
    arma::vec nv = nr.nv();

    ofs << i << '\n';
    for (arma::uword j = 0; j < tot; j++) {
      if (j == i && j < (tot - 1))
      {
        ofs << 0 << '\t' << nv(j) << '\t';
      }
      else if (i == j && j == (tot - 1))
      {
        ofs << 0 << '\t' << nv(j) << '\n';
      }
      else if (i > j && j == (tot - 1))
      {
        ofs << nv(j) << '\t' << 0 << '\n';
      }
      else
      {
        ofs << nv(j) <<
            (j == (tot - 1) ? '\n' : '\t');
      }
    }
    time = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    log << "Time for gene " << genes[i] << "\t" << time << '\n';
    log.send(LOG_INFO);
  }
  ofs.close();
}
