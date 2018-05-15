// Seidr
#include <common.h>
#include <elnet-fun.h>
#include <mpims.h>
#include <glmnetx.h>
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
#include <cmath>

namespace fs = boost::filesystem;

std::random_device rd;
std::mt19937 gen(rd());

void print_null(const char *s) {}

class seidr_mpi_elnet : public seidr_mpi {
public:
  void entrypoint();
  void finalize();
  seidr_mpi_elnet() : seidr_mpi() {}
  seidr_mpi_elnet(unsigned long bs) : seidr_mpi(bs) {}
  void set_alpha(double x) {_alpha = x;}
  void set_flmin(double x) {_flmin = x;}
  void set_nlam(arma::uword x) {_nlam = x;}
  void set_min_sample_size(arma::uword x) {_min_sample_size = x;}
  void set_max_sample_size(arma::uword x) {_max_sample_size = x;}
  void set_predictor_sample_size_min(arma::uword x) {_predictor_sample_size_min = x;}
  void set_predictor_sample_size_max(arma::uword x) {_predictor_sample_size_max = x;}
  void set_ensemble_size(arma::uword x) {_ensemble_size = x;}
  void set_targeted(bool x) { _targeted = x; }
private:
  arma::uword _min_sample_size = 0;
  arma::uword _max_sample_size = 0;
  arma::uword _predictor_sample_size_min = 0;
  arma::uword _predictor_sample_size_max = 0;
  arma::uword _ensemble_size = 0;
  double _alpha = 0;
  double _flmin = 0;
  arma::uword _nlam = 0;
  bool _targeted = false;
};

void seidr_mpi_elnet::entrypoint()
{
  seidr_mpi_logger log;
  std::string tmpfile = _outfilebase + "/.seidr_tmp_elnet/MPIthread_" +
                        std::to_string(_id) + "_" + std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);

  std::vector<arma::uword> uvec;
  for (auto i : _indices)
    uvec.push_back(i);
  el_ensemble(_data, _genes, uvec, tmpfile, _id, _min_sample_size,
              _max_sample_size, _predictor_sample_size_min,
              _predictor_sample_size_max, _ensemble_size, _alpha, _flmin,
              _nlam);
  announce_ready();
}

void seidr_mpi_elnet::finalize()
{
  merge_files(_outfile, _outfilebase, ".seidr_tmp_elnet",
              _targeted, _id, _genes);
  check_logs();
}

void el_ensemble(arma::mat geneMatrix, std::vector<std::string> genes,
                 std::vector<arma::uword> uvec, std::string outfile,
                 int thread_id, arma::uword min_sample_size,
                 arma::uword max_sample_size,
                 arma::uword predictor_sample_size_min,
                 arma::uword predictor_sample_size_max,
                 arma::uword ensemble_size, double alpha, double flmin,
                 arma::uword nlam) {

  seidr_mpi_logger log;

  std::ofstream ofs(outfile.c_str(), std::ios::out);
  if (! ofs)
    throw std::runtime_error("Could not open temp file: " + outfile);

  arma::vec ret(geneMatrix.n_cols);

  std::uniform_int_distribution<> sample_gen(min_sample_size,
      max_sample_size);
  std::uniform_int_distribution<> predictor_gen(predictor_sample_size_min,
      predictor_sample_size_max);

  seidr_score_t fensemble_size = 0;

  try
  {
    fensemble_size = boost::numeric_cast<seidr_score_t>(ensemble_size);
  }
  catch (boost::numeric::bad_numeric_cast& e)
  {
    throw std::runtime_error(e.what());
  }


  arma::uvec samples(geneMatrix.n_rows);
  for (arma::uword i = 0; i < geneMatrix.n_rows; i++)
  {
    samples(i) = i;
  }

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
    // Ensure all genes have the same initial seed such that
    // sample sizes are identical
    std::seed_seq sseq{3, 1, 4, 1, 5, 9, 2, 6, 5};
    gen.seed(sseq);
    arma::uword reshuf = 0;
    for (arma::uword boot = 0; boot < ensemble_size;)
    {
      try
      {
        // Randomly select a number of predictor genes in range
        // predictor_sample_size_min to predictor_sample_size_max
        pred = arma::shuffle(pred);
        arma::uword pred_size = predictor_gen(gen);
        arma::uvec pred_sub(pred_size);
        for (arma::uword i = 0; i < pred_size; i++)
        {
          pred_sub(i) = pred(i);
        }

        // Randomly select a number of samples in range
        // min_sample_size to max_sample_size
        samples = arma::shuffle(samples);
        arma::uword sample_size = sample_gen(gen);
        arma::uvec sample_sub(sample_size);
        for (arma::uword i = 0; i < sample_size; i++)
        {
          sample_sub(i) = samples(i);
        }

        log << "Predictors: " << pred_size << ", "
            << "Samples: " << sample_size << ", "
            << "Iteration: " << boot << '\n';
        log.send(LOG_DEBUG);
        // Create subset matrix
        arma::mat pred_mat = geneMatrix.submat(sample_sub, pred_sub);

        // Get target expression and make sure we have the same samples
        // as in the predictor expression
        arma::vec target_exp = geneMatrix.col(target);
        target_exp = target_exp.elem(sample_sub);

        glm elnet(pred_mat, target_exp, nlam, flmin, alpha);

        double nr = pred_mat.n_rows;
        arma::uword k = round(sqrt(nr));

        if (k > 10) k = 10;

        glm_cv_t glm_cv = elnet.k_fold_cv(k);

        arma::vec lambda{glm_cv.lambda_min};

        glm best_elnet(pred_mat, target_exp, nlam, 2, alpha, lambda);
        best_elnet.calculate_beta();

        // Calculate absolute value of weights and rank
        arma::mat weights = best_elnet.beta;
        weights = arma::abs(weights);
        arma::uvec sort_index = arma::sort_index(weights, "descend");

        // Check which weights were zero
        std::vector<arma::uword> zeros;
        for (arma::uword i = 0; i < weights.n_elem; i++)
        {
          if (almost_equal(weights(i), 0))
          {
            zeros.push_back(i);
          }
        }
        arma::uvec zero_ind = arma::conv_to<arma::uvec>::from(zeros);

        //weights.print();
        // Zero out weights vector and set top 20th values to 1
        weights.zeros();
        arma::uword rt = pred_mat.n_cols / 20;
        for (arma::uword i = 0; i <= rt; i++)
        {
          weights(sort_index(i)) = 1;
        }

        // Reset those weights which were zero before ranking
        weights.elem(zero_ind).zeros();

        // Update return vector with new model data
        ret.elem(pred_sub) += weights;
        boot++;
      }
      catch (std::exception& e)
      {
        reshuf++;
        if (reshuf > boot / 3)
        {
          log << "Gene " << genes[target] << " was too low variance " <<
              "and ancountered too many reshuffles. All output will be 0.\n";
          log.send(LOG_WARN);
          ret.zeros();
          break;
        }
      }
    }

    ofs << target << '\n';
    for (seidr_uword_t i = 0; i < ret.size(); i++)
    {
      ofs << ret[i] / fensemble_size
          << (i == ret.size() - 1 ? '\n' : '\t');
    }
  }

}

void el_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
             std::string outfile, arma::uword min_sample_size,
             arma::uword max_sample_size,
             arma::uword predictor_sample_size_min,
             arma::uword predictor_sample_size_max,
             arma::uword ensemble_size, double alpha, double flmin,
             arma::uword nlam) {

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<unsigned long> uvec;
  for (unsigned long i = 0; i < GM.n_cols; i++)
    uvec.push_back(i);
  seidr_mpi_elnet mpi(bs);

  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(uvec);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_min_sample_size(min_sample_size);
  mpi.set_max_sample_size(max_sample_size);
  mpi.set_predictor_sample_size_min(predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(predictor_sample_size_max);
  mpi.set_ensemble_size(ensemble_size);
  mpi.set_alpha(alpha);
  mpi.set_flmin(flmin);
  mpi.set_nlam(nlam);

  mpi.exec();
}

void el_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                std::vector<std::string> targets, std::string outfile,
                arma::uword min_sample_size,
                arma::uword max_sample_size,
                arma::uword predictor_sample_size_min,
                arma::uword predictor_sample_size_max,
                arma::uword ensemble_size,
                double alpha, double flmin,
                arma::uword nlam) {

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

  seidr_mpi_elnet mpi(bs);
  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(positions);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_min_sample_size(min_sample_size);
  mpi.set_max_sample_size(max_sample_size);
  mpi.set_predictor_sample_size_min(predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(predictor_sample_size_max);
  mpi.set_ensemble_size(ensemble_size);
  mpi.set_alpha(alpha);
  mpi.set_flmin(flmin);
  mpi.set_nlam(nlam);
  mpi.set_targeted(true);

  mpi.exec();
}
