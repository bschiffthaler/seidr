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

// Seidr
#include <BSlogger.h>
#include <common.h>
#include <elnet-fun.h>
#include <fs.h>
#include <glmnetx.h>
#include <mpiomp.h>

#include <omp.h>
// External
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>

#include <armadillo>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

using boost::numeric_cast;

void print_null(const char *s) {}

class seidr_mpi_elnet : public seidr_mpi_omp {
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_alpha(double x) {_alpha = x;}
  void set_flmin(double x) {_flmin = x;}
  void set_nlam(arma::uword x) {_nlam = x;}
  void set_min_sample_size(arma::uword x) {_min_sample_size = x;}
  void set_max_sample_size(arma::uword x) {_max_sample_size = x;}
  void set_predictor_sample_size_min(arma::uword x) {
    _predictor_sample_size_min = x;
  }
  void set_predictor_sample_size_max(arma::uword x) {
    _predictor_sample_size_max = x;
  }
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
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  while (! _my_indices.empty())
  {
    std::vector<arma::uword> uvec;
    for (auto i : _my_indices)
    {
      uvec.push_back(i);
    }

    el_ensemble(_data, _genes, uvec, _tempdir, _min_sample_size,
                _max_sample_size, _predictor_sample_size_min,
                _predictor_sample_size_max, _ensemble_size, _alpha, _flmin,
                _nlam, this);
    get_more_work();
  }
  log << "No more work. Waiting for other tasks to finish...\n";
  log.send(LOG_INFO);
}

void seidr_mpi_elnet::finalize()
{
  merge_files(_outfile, _tempdir, _targeted, _id, _genes);
}

void el_ensemble(const arma::mat& geneMatrix,
                 const std::vector<std::string>& genes,
                 const std::vector<arma::uword>& uvec,
                 const std::string& tmpdir,
                 const arma::uword min_sample_size,
                 const arma::uword max_sample_size,
                 const arma::uword predictor_sample_size_min,
                 const arma::uword predictor_sample_size_max,
                 const arma::uword ensemble_size,
                 const double alpha,
                 const double flmin,
                 const arma::uword nlam,
                 seidr_mpi_elnet* self) {

  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<> sample_gen(min_sample_size,
      max_sample_size);
  std::uniform_int_distribution<> predictor_gen(predictor_sample_size_min,
      predictor_sample_size_max);

  seidr_score_t fensemble_size = 0;


  std::string tmpfile = tempfile(tmpdir);
  std::ofstream ofs(tmpfile.c_str(), std::ios::out);


  fensemble_size = boost::numeric_cast<seidr_score_t>(ensemble_size);

  arma::uvec samples(geneMatrix.n_rows);
  for (arma::uword i = 0; i < geneMatrix.n_rows; i++)
  {
    samples(i) = i;
  }

  #pragma omp parallel for
  for (uint64_t i = 0; i < uvec.size(); i++)  // NOLINT 
  {
    arma::vec ret(geneMatrix.n_cols);
    seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
    const arma::uword& target = uvec[i];
    #pragma omp critical
    {
      log << "Started gene: " << genes[target] << ".\n";
      log.send(LOG_INFO);
      while (self->check_logs(LOG_NAME"@" + mpi_get_host())); // NOLINT
    }
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

        // Create subset matrix
        arma::mat pred_mat = geneMatrix.submat(sample_sub, pred_sub);

        // Get target expression and make sure we have the same samples
        // as in the predictor expression
        arma::vec target_exp = geneMatrix.col(target);
        target_exp = target_exp.elem(sample_sub);

        glm elnet(pred_mat, target_exp, nlam, flmin, alpha);

        double nr = pred_mat.n_rows;
        auto k = numeric_cast<arma::uword>(round(sqrt(nr)));

        if (k > 10)
        {
          k = 10;
        }

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
        bool should_break = false;
        #pragma omp critical
        {
          if (reshuf > boot / 3)
          {
            log << "Gene " << genes[target] << " was too low variance " <<
                "and ancountered too many reshuffles. All output will be 0.\n";
            log.send(LOG_WARN);
            ret.zeros();
            should_break = true;
          }
        }
        if (should_break)
        {
          break;
        }
      }
    }
    #pragma omp critical
    {
      ofs << target << '\n';
      for (seidr_uword_t i = 0; i < ret.size(); i++)
      {
        ofs << ret[i] / fensemble_size
            << (i == ret.size() - 1 ? '\n' : '\t');
      }
    }
  }
  ofs.close();
}

void el_full(const arma::mat& GM,
             const std::vector<std::string>& genes,
             const seidr_elnet_param_t& param) {
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  std::vector<uint64_t> uvec;
  for (uint64_t i = 0; i < GM.n_cols; i++)
  {
    uvec.push_back(i);
  }

  seidr_mpi_elnet mpi(param.bs, GM, uvec  , genes, param.tempdir,
                      param.outfile);
  mpi.set_min_sample_size(param.min_sample_size);
  mpi.set_max_sample_size(param.max_sample_size);
  mpi.set_predictor_sample_size_min(param.predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(param.predictor_sample_size_max);
  mpi.set_ensemble_size(param.ensemble_size);
  mpi.set_alpha(param.alpha);
  mpi.set_flmin(param.flmin);
  mpi.set_nlam(param.nlam);

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
    }
  }

  MPI_Finalize();
}

void el_partial(const arma::mat& GM,
                const std::vector<std::string>& genes,
                const std::vector<std::string>& targets,
                const seidr_elnet_param_t& param) {
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  std::vector<uint64_t> positions;
  for (uint64_t i = 0; i < targets.size(); i++) {
    uint64_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size())
    {
      #pragma omp critical
      {
        log << "Gene " << targets[i]
        << " was not found in the expression set "
        << "and will therefore not be considered."
        << " Please check that your expression set and "
        << "its column names (gene file) contain an entry for "
        << targets[i] << ".\n";
        log.send(LOG_WARN);
      }
    }
    else
    {
      positions.push_back(pos);
    }
  }

  seidr_mpi_elnet mpi(param.bs, GM, positions, genes, param.tempdir,
                      param.outfile);
  mpi.set_min_sample_size(param.min_sample_size);
  mpi.set_max_sample_size(param.max_sample_size);
  mpi.set_predictor_sample_size_min(param.predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(param.predictor_sample_size_max);
  mpi.set_ensemble_size(param.ensemble_size);
  mpi.set_alpha(param.alpha);
  mpi.set_flmin(param.flmin);
  mpi.set_nlam(param.nlam);
  mpi.set_targeted(true);

  mpi.entrypoint();

  MPI_Barrier(MPI_COMM_WORLD); // NOLINT
  mpi.remove_queue_file();
  if (mpi.rank() == 0)
  {
    log << "Finalizing...\n";
    log.send(LOG_INFO);
    mpi.finalize();
  }

  MPI_Finalize();
}
