// Seidr
#include <common.h>
#include <plsnet-fun.h>
#include <mpims.h>
// External
#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <map>
#include <armadillo>
#include <ctime>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace fs = boost::filesystem;

typedef std::pair<std::streampos, std::string> file_index_t;
typedef std::map<seidr_uword_t, file_index_t> result_file_map_t;

std::random_device rd;
std::mt19937 gen(rd());

class seidr_mpi_plsnet : public seidr_mpi {
public:
  void entrypoint();
  void finalize();
  seidr_mpi_plsnet() : seidr_mpi() {}
  seidr_mpi_plsnet(unsigned long bs) : seidr_mpi(bs) {}
  void set_predictor_sample_size(arma::uword x) {_predictor_sample_size = x;}
  void set_ensemble_size(arma::uword x) {_ensemble_size = x;}
  void set_ncomp(arma::uword x) {_ncomp = x;}
  void set_targeted(bool x) {_targeted = x;}
private:
  arma::uword _predictor_sample_size;
  arma::uword _ensemble_size;
  arma::uword _ncomp;
  bool _targeted;
};

void seidr_mpi_plsnet::entrypoint()
{
  seidr_mpi_logger log;
  std::string tmpfile = _outfilebase + "/.seidr_tmp_plsnet/MPIthread_" +
                        std::to_string(_id) + "_" + std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);

  std::vector<arma::uword> uvec;
  for (auto i : _indices)
    uvec.push_back(i);
  plsnet(_data, _genes, uvec, tmpfile, _id, _predictor_sample_size,
         _ensemble_size, _ncomp);
  announce_ready();
}

void seidr_mpi_plsnet::finalize()
{
  result_file_map_t rmap;
  if (_id == 0)
  {
    seidr_mpi_logger log;
    log << "Merging tmp files and cleaning up.\n";
    log.send(LOG_INFO);
    std::vector<fs::path> files;
    std::string tmpdir = _outfilebase + "/.seidr_tmp_plsnet/";
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

    arma::mat VIM(rmap.size(), _data.n_cols, arma::fill::zeros);
    arma::uword i = 0;
    auto it = rmap.begin();
    auto gene_index = it->second.first;
    auto file_path = it->second.second;
    std::ifstream ifs(file_path);
    for (; it != rmap.end(); it++)
    {
      ifs.seekg(gene_index);
      std::string l;
      std::getline(ifs, l);
      std::stringstream ss(l);
      std::string field;
      arma::uword j = 0;
      while (std::getline(ss, field, '\t'))
      {
        double value = std::stod(field);
        VIM(i, j) = value;
        j++;
      }
      i++;
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

    // Refine matrix
    arma::rowvec theta = arma::var(VIM, 0, 0);
    for (arma::uword col = 0; col < VIM.n_cols; col++)
      VIM.col(col) *= theta(col);

    if (_targeted)
    {
      arma::uword ctr = 0;
      for (auto tar : rmap)
      {
        arma::uword row = tar.first;
        for(arma::uword col = 0; col < VIM.n_cols; col++)
        {
          if(row == col) continue;
          ofs << _genes[row] << '\t' << _genes[col] << '\t' << VIM(ctr, col) << '\n';
        }
        ctr++;
      }
    }
    else
    {
      for (arma::uword row = 0; row < VIM.n_rows; row++)
      {
        for (arma::uword col = 0; col < VIM.n_cols; col++)
        {
          ofs << VIM(row, col)
              << (col == VIM.n_cols - 1 ? '\n' : '\t');
        }
      }
    }

    fs::remove_all(p_tmp);
    log << "Waiting up to 10 seconds for queued logs...\n";
    log.send(LOG_INFO);
    double now = MPI_Wtime();
    // Push all pending logs
    while (MPI_Wtime() - now < 10)
    {
      check_logs();
    }
  }
}


void plsnet(arma::mat geneMatrix, std::vector<std::string> genes,
            std::vector<arma::uword> uvec, std::string outfile,
            int thread_id,
            arma::uword predictor_sample_size,
            arma::uword ensemble_size, arma::uword ncomp)
{

  seidr_mpi_logger log;
  std::ofstream ofs(outfile.c_str(), std::ios::out);
  if (! ofs)
    throw std::runtime_error("Could not open temp file: " + outfile);

  // Create random generators for samples and genes
  std::uniform_int_distribution<> sample_gen(0, geneMatrix.n_rows - 1);

  arma::vec ret(geneMatrix.n_cols, arma::fill::zeros);

  for (auto& target : uvec)
  {
    log << "Started gene: " << genes[target] << ".\n";
    log.send(LOG_INFO);
    arma::uvec pred(geneMatrix.n_cols - 1);
    arma::uword j = 0;
    for (arma::uword i = 0; i < geneMatrix.n_cols; i++)
    {
      if (i != target)
      {
        pred(j) = i;
        j++;
      }
    }

    ret.zeros();
    std::seed_seq seeds{3, 1, 4, 1, 5, 9, 2, 6, 5};
    gen.seed(seeds);
    arma::arma_rng::set_seed(314159265);

    for (arma::uword i = 0; i < ensemble_size; i++)
    {
      // Generate vector of random samples with replacement
      arma::uvec samples(geneMatrix.n_rows);
      for (arma::uword j = 0; j < geneMatrix.n_rows; j++)
      {
        samples(j) = sample_gen(gen);
      }

      // Generate vector of random predictors without replacement
      arma::uvec pred_sub(predictor_sample_size);
      pred = arma::shuffle(pred);
      for (arma::uword j = 0; j < predictor_sample_size; j++)
      {
        pred_sub(j) = pred(j);
      }

      // Subset expression matrix
      arma::mat X = geneMatrix.submat(samples, pred_sub);
      arma::vec Y = geneMatrix.col(target);
      Y = Y(samples);
      // Get variable importance
      arma::vec vim = vip(X, Y, ncomp);
      ret(pred_sub) += vim;
    }
    ofs << target << '\n';
    for (seidr_uword_t i = 0; i < ret.size(); i++)
    {
      ofs << ret[i] << (i == ret.size() - 1 ? '\n' : '\t');
    }
  }
}

void plsnet_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
                 std::string outfile,
                 arma::uword predictor_sample_size,
                 arma::uword ensemble_size, arma::uword ncomp) {

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<unsigned long> uvec;
  for (unsigned long i = 0; i < GM.n_cols; i++)
    uvec.push_back(i);
  seidr_mpi_plsnet mpi(bs);

  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(uvec);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_predictor_sample_size(predictor_sample_size);
  mpi.set_ensemble_size(ensemble_size);
  mpi.set_ncomp(ncomp);

  mpi.exec();
}

void plsnet_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                    std::vector<std::string> targets, std::string outfile,
                    arma::uword predictor_sample_size,
                    arma::uword ensemble_size, arma::uword ncomp) {

  seidr_mpi_logger log;

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<size_t> positions;
  for (size_t i = 0; i < targets.size(); i++) {
    size_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size()) {
      log << "Gene " << targets[i]
          << " was not found in the expression set "
          << "and will therefore not be considered."
          << " Please check that your expression set and "
          << "its column names (gene file) contain an entry for "
          << targets[i] << ".\n";
      log.log(LOG_WARN);
    } else {
      positions.push_back(pos);
    }
  }

  seidr_mpi_plsnet mpi(bs);
  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(positions);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_predictor_sample_size(predictor_sample_size);
  mpi.set_ensemble_size(ensemble_size);
  mpi.set_ncomp(ncomp);
  mpi.set_targeted(true);

  mpi.exec();
}

arma::vec vip(arma::mat& X, arma::vec& Y, arma::uword ncomp)
{
  plsreg_t plsr = plsreg(X, Y, ncomp);
  arma::mat wsq = arma::pow(plsr.w, 2);
  wsq = wsq.t();
  double xnr = X.n_cols;
  arma::mat xret = xnr * plsr.p.row(1) * wsq / arma::accu(plsr.p.row(1));
  xret = xret.t();
  return xret.col(0);
}

plsreg_t plsreg(arma::mat& X, arma::vec& Y, arma::uword ncomp)
{
  plsreg_t ret;

  arma::mat X0 = X;
  X0.each_col([](arma::vec & v) { v -= arma::mean(v); });

  arma::vec Y0 = Y - arma::mean(Y);

  simpls_t spls = simpls(X0, Y0, ncomp);

  arma::mat pctvar(2, ncomp, arma::fill::zeros);

  pctvar.row(0) = arma::sum(arma::pow(spls.x, 2), 0) /
                  arma::accu(arma::sum(arma::pow(X0, 2), 0));
  pctvar.row(1) = arma::sum(arma::pow(spls.y, 2), 0) /
                  arma::accu(arma::sum(arma::pow(Y0, 2), 0));

  ret.p = pctvar;
  ret.w = spls.w;

  return ret;
}

simpls_t simpls(arma::mat& X, arma::vec& Y, arma::uword ncomp)
{
  simpls_t ret;
  arma::mat Cov = X.t() * Y;

  arma::mat Xloadings(X.n_cols, ncomp, arma::fill::zeros);
  arma::rowvec Yloadings(ncomp, arma::fill::zeros);
  arma::mat Weights = Xloadings;
  arma::mat W = Xloadings;
  //cov.print();
  for (arma::uword i = 0; i < ncomp; i++)
  {
    arma::mat U, V;
    arma::vec s;
    svd_econ(U, s, V, Cov);
    U = U.col(0);
    s = s(0);
    V = V.col(0);
    arma::vec ti = X * U;
    double normti = norm(ti);
    ti /= normti;
    Xloadings.col(i) = X.t() * ti;
    double qi = s(0) * V(0, 0) / normti;
    Yloadings(i) = qi;
    U /= normti;
    Weights.col(i) = U;
    arma::vec vi = Xloadings.col(i);
    for (arma::uword r = 0; r < 2; r++)
    {
      for (arma::uword j  = 0; j < i; j++)
      {
        arma::vec vj = W.col(j);
        arma::vec tmp = (vi.t() * vj);
        vi -= vj * tmp(0);
      };
    }
    vi /= norm(vi);
    W.col(i) = vi;
    Cov = Cov - vi * (vi.t() * Cov);
    arma::mat Vi = W.cols(0, i);
    Cov = Cov - Vi * (Vi.t() * Cov);
  }
  ret.x = Xloadings;
  ret.y = Yloadings;
  ret.w = Weights;
  return ret;
}
