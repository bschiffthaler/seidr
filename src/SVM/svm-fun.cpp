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
#include <common.h>
#include <svm-fun.h>
#ifdef SEIDR_WITH_MPI
  #include <mpiomp.h>
#else
  #include <mpi_dummy.h>
#endif
// External
#include <iostream>
#include <random>
#include <string>
#include <svm.h>
#include <fstream>
#include <map>
#include <armadillo>
#include <ctime>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace fs = boost::filesystem;

std::random_device rd;
std::mt19937 gen(rd());

void print_null(const char *s) {}

class seidr_mpi_svm : public seidr_mpi_omp {
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_param(svm_parameter& p) {_param = p;}
  void set_min_sample_size(arma::uword x) {_min_sample_size = x;}
  void set_max_sample_size(arma::uword x) {_max_sample_size = x;}
  void set_predictor_sample_size_min(arma::uword x) {_predictor_sample_size_min = x;}
  void set_predictor_sample_size_max(arma::uword x) {_predictor_sample_size_max = x;}
  void set_ensemble_size(arma::uword x) {_ensemble_size = x;}
  void set_targeted(bool x) {_targeted = x;}
private:
  svm_parameter _param;
  arma::uword _min_sample_size = 0;
  arma::uword _max_sample_size = 0;
  arma::uword _predictor_sample_size_min = 0;
  arma::uword _predictor_sample_size_max = 0;
  arma::uword _ensemble_size = 0;
  bool _targeted = false;
};

void seidr_mpi_svm::entrypoint()
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  while (! _my_indices.empty())
  {
    std::vector<arma::uword> uvec;
    for (auto i : _my_indices)
    {
      uvec.push_back(i);
    }
    svm(_data, _genes, uvec, _tempdir, _param, _min_sample_size,
        _max_sample_size, _predictor_sample_size_min,
        _predictor_sample_size_max, _ensemble_size, this);
    get_more_work();
  }
  log << "No more work. Waiting for other tasks to finish...\n";
  log.send(LOG_INFO);
}

void seidr_mpi_svm::finalize()
{
  merge_files(_outfile, _tempdir,
              _targeted, _id, _genes);
}

void svm(const arma::mat& geneMatrix,
         const std::vector<std::string>& genes,
         const std::vector<arma::uword>& uvec,
         const std::string& tmpdir,
         svm_parameter& param,
         const arma::uword& min_sample_size,
         const arma::uword& max_sample_size,
         const arma::uword& predictor_sample_size_min,
         const arma::uword& predictor_sample_size_max,
         const arma::uword& ensemble_size,
         seidr_mpi_svm * self) {

  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  std::string tmpfile = tempfile(tmpdir);

  std::ofstream ofs(tmpfile.c_str(), std::ios::out);
  svm_set_print_string_function(print_null);
  if (! ofs)
  {
    throw std::runtime_error("Could not open temp file: " + tmpfile);
  }

  std::uniform_int_distribution<>
  sample_gen(min_sample_size, max_sample_size);
  std::uniform_int_distribution<>
  predictor_gen(predictor_sample_size_min, predictor_sample_size_max);

  seidr_score_t fensemble_size = 0;

  fensemble_size = boost::numeric_cast<seidr_score_t>(ensemble_size);

  arma::uvec samples(geneMatrix.n_rows);
  for (arma::uword i = 0; i < geneMatrix.n_rows; i++)
  {
    samples(i) = i;
  }

  #pragma omp parallel for schedule(dynamic)
  for (uint64_t i = 0; i < uvec.size(); i++)
  {
    arma::vec ret(geneMatrix.n_cols);
    auto& target = uvec[i];
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
    arma::arma_rng::set_seed(314159265);
    for (arma::uword boot = 0; boot < ensemble_size; boot++)
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

      // log << "Predictors: " << pred_size << ", "
      //     << "Samples: " << sample_size << ", "
      //     << "Iteration: " << boot << '\n';
      // log.send(LOG_DEBUG);
      // Create subset matrix
      arma::mat pred_mat = geneMatrix.submat(sample_sub, pred_sub);

      // Set up SVM problem struct sparse node notation
      svm_node ** nodes = new svm_node*[pred_mat.n_rows];

      for (arma::uword i = 0; i < pred_mat.n_rows; i++)
      {
        nodes[i] = new svm_node[pred_mat.n_cols + 1];
      }

      for (arma::uword i = 0; i < pred_mat.n_rows; i++)
      {
        arma::uword t = 0;
        for (arma::uword j = 0; j < pred_mat.n_cols; j++)
        {
          int x = 0;
          // Safely cast uword to int
          try
          {
            x = boost::numeric_cast<int>(j);
          }
          catch (boost::numeric::bad_numeric_cast& e)
          {
            throw std::runtime_error(e.what());
          }
          double y = pred_mat(i, j);
          // Do not add zeros to the problem
          if (! almost_equal(y, 0))
          {
            nodes[i][t] = svm_node{x, y};
            t++;
          }
        }
        // The end of the nodes is signified by a node with index -1
        nodes[i][t] = svm_node{ -1, 0};
      }

      // Get target expression and make sure we have the same samples
      // as in the predictor expression
      arma::vec target_exp = geneMatrix.col(target);
      target_exp = target_exp.elem(sample_sub);

      int nrows = 0;

      try
      {
        nrows = boost::numeric_cast<int>(pred_mat.n_rows);
      }
      catch (boost::numeric::bad_numeric_cast& e)
      {
        throw std::runtime_error(e.what());
      }

      svm_problem prob{nrows, target_exp.memptr(), nodes};

      // Check parameters of the model
      const char * error_msg = svm_check_parameter(&prob, &param);
      if (error_msg)
      {
        std::string msg(error_msg);
        throw std::runtime_error(msg);
      }

      svm_model * model = svm_train(&prob, &param);

      // Extract support vector coefficients and indices
      arma::mat coefs(1, model->l);
      for (int i = 0; i < model->l; i++)
      {
        coefs(0, i) = model->sv_coef[0][i];
      }

      arma::uvec sv_ind(model->l);
      int nr_sv = svm_get_nr_sv(model);
      int * sv_indices = new int[nr_sv];
      svm_get_sv_indices(model, sv_indices);
      for (int i = 0; i < nr_sv; i++)
      {
        sv_ind[i] = sv_indices[i] - 1;
      }

      // Get support vectors from the predictor matrix
      arma::mat svs = pred_mat.rows(sv_ind);

      // Calculate absolute value of weights and rank
      arma::mat weights = coefs * svs;
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

      // Free memory of SVM problem and model
      svm_free_and_destroy_model(&model);
      delete[] sv_indices;
      for (arma::uword i = 0; i < pred_mat.n_rows; i++)
      {
        delete[] nodes[i];
      }
      delete[] nodes;
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

void svm_full(const arma::mat& GM,
              const std::vector<std::string>& genes, 
              seidr_svm_param_t& param) {
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> uvec;
  for (uint64_t i = 0; i < GM.n_cols; i++)
  {
    if (param.resuming)
    {
      if (! in_sorted_range<uint64_t>(i, param.good_idx))
      {
        uvec.push_back(i);
      }
    }
    else
    {
      uvec.push_back(i);
    }
  }

  seidr_mpi_svm mpi(param.bs, GM, uvec, genes, param.tempdir,
                       param.outfile);
  mpi.set_param(param.svparam);
  mpi.set_min_sample_size(param.min_sample_size);
  mpi.set_max_sample_size(param.max_sample_size);
  mpi.set_predictor_sample_size_min(param.predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(param.predictor_sample_size_max);
  mpi.set_ensemble_size(param.ensemble_size);

  mpi.entrypoint();

  SEIDR_MPI_BARRIER(); // NOLINT
  svm_destroy_param(&param.svparam);
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

  SEIDR_MPI_FINALIZE();
}

void svm_partial(const arma::mat& GM,
                 const std::vector<std::string>& genes,
                 const std::vector<std::string>& targets,
                 seidr_svm_param_t& param) {

  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> positions;
  for (uint64_t i = 0; i < targets.size(); i++) {
    uint64_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size()) 
    {
      log << "Gene " << targets[i]
          << " was not found in the expression set "
          << "and will therefore not be considered."
          << " Please check that your expression set and "
          << "its column names (gene file) contain an entry for "
          << targets[i] << ".\n";
      log.log(LOG_WARN);
    } 
    else 
    {
      if (param.resuming)
      {
        if (! in_sorted_range<uint64_t>(pos, param.good_idx))
        {
          positions.push_back(pos);
        }
      }
      else
      {
        positions.push_back(pos);
      }
    }
  }

  seidr_mpi_svm mpi(param.bs, GM, positions, genes, param.tempdir,
                       param.outfile);
  mpi.set_param(param.svparam);
  mpi.set_min_sample_size(param.min_sample_size);
  mpi.set_max_sample_size(param.max_sample_size);
  mpi.set_predictor_sample_size_min(param.predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(param.predictor_sample_size_max);
  mpi.set_ensemble_size(param.ensemble_size);
  mpi.set_targeted(true);
  mpi.entrypoint();

  SEIDR_MPI_BARRIER(); // NOLINT
  svm_destroy_param(&param.svparam);
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

  SEIDR_MPI_FINALIZE();
}
