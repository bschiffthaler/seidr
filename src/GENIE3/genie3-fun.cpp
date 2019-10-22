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
#include <genie3-fun.h>
#include <mpiomp.h>
// External
#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <armadillo>
#include <ctime>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <ranger/DataDouble.h>
#include <ranger/utility.h>
#include <ranger/ForestRegression.h>
#include <ranger/globals.h>
#include <map>
#include <sstream>
#include <memory>

namespace fs = boost::filesystem;

class seidr_mpi_genie3 : public seidr_mpi_omp {
public:
  using seidr_mpi_omp::seidr_mpi_omp;
  void entrypoint();
  void finalize();
  void set_ntree(uint64_t x) {_ntree = x;}
  void set_mtry(uint64_t x) {_mtry = x;}
  void set_min_node_size(uint64_t x) {_min_node_size = x;}
  void set_alpha(double x) {_alpha = x;}
  void set_minprop(double x) {_minprop = x;}
  void set_targetd(bool x) {_targeted = x;}
private:
  uint64_t _ntree = 0;
  uint64_t _mtry = 0;
  uint64_t _min_node_size = 0;
  double _alpha = 0;
  double _minprop = 0;
  bool _targeted = 0;
};

void seidr_mpi_genie3::entrypoint()
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  while (! _my_indices.empty())
  {
    std::vector<arma::uword> uvec;
    for (auto i : _my_indices)
    {
      uvec.push_back(i);
    }
    genie3(_data, _genes, uvec, _tempdir, _ntree, _mtry,
           _min_node_size, _alpha, _minprop);
    get_more_work();
  }
  log << "No more work. Waiting for other tasks to finish...\n";
  log.send(LOG_INFO);
}

void seidr_mpi_genie3::finalize()
{
  merge_files(_outfile, _tempdir, _targeted, _id, _genes);
}

class SeidrForestData : public ranger::Data {
public:
  SeidrForestData() = default;

  SeidrForestData(const SeidrForestData&) = delete;
  SeidrForestData& operator=(const SeidrForestData&) = delete;

  virtual ~SeidrForestData() override = default;

  double get_x(size_t row, size_t col) const override {
    // Use permuted data for corrected impurity importance
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    if (col < num_cols_no_snp) {
      return x[col * num_rows + row];
    } else {
      return getSnp(row, col, col_permuted);
    }
  }

  double get_y(size_t row, size_t col) const override {
    return y[col * num_rows + row];
  }

  void reserveMemory(size_t y_cols) override {
    x.resize(num_cols * num_rows);
    y.resize(y_cols * num_rows);
  }

  void set_x(size_t col, size_t row, double value, bool& error) override {
    x[col * num_rows + row] = value;
  }

  void set_y(size_t col, size_t row, double value, bool& error) override {
    y[col * num_rows + row] = value;
  }


private:
  std::vector<double> x;
  std::vector<double> y;
};

class SeidrForest : public ranger::ForestRegression {
public:
  SeidrForest(const std::vector<std::string>& dependent_variable_name,
              const arma::mat& gm_x,
              const arma::mat& gm_y,
              const std::vector<std::string>& var,
              const size_t ntree,
              const size_t mtry,
              const size_t min_node_size,
              const double alpha,
              const double minprop,
              std::ostream& ostream);
  std::unique_ptr<ranger::Data> load(const std::vector<std::string>& dependent_variable_name,
                                     const arma::mat& gm_x,
                                     const arma::mat& gm_y,
                                     const std::vector<std::string>& var);
};

std::unique_ptr<ranger::Data> SeidrForest::load(const std::vector<std::string>& dependent_variable_name,
    const arma::mat& gm_x,
    const arma::mat& gm_y,
    const std::vector<std::string>& var)
{
  std::unique_ptr<ranger::Data> internal_data_ptr(new SeidrForestData());
  internal_data_ptr->load_seidr(dependent_variable_name, gm_x, gm_y, var);
  return internal_data_ptr;
}

SeidrForest::SeidrForest(const std::vector<std::string>& dependent_variable_name,
                         const arma::mat& gm_x,
                         const arma::mat& gm_y,
                         const std::vector<std::string>& var,
                         const size_t ntree,
                         const size_t mtry,
                         const size_t min_node_size,
                         const double alpha,
                         const double minprop,
                         std::ostream& ostream) {

  std::vector<std::string> catvars;

  ranger::MemoryMode M = ranger::MEM_DOUBLE;
  ranger::ImportanceMode I = ranger::IMP_GINI;
  ranger::SplitRule S = ranger::LOGRANK;
  ranger::PredictionType P = ranger::RESPONSE;

  std::vector<double> sample_fraction_vector = { 1 };

  // std::vector<double> internal_data;
  // for (size_t i = 0; i < gm.n_cols; i++)
  // {
  //   for (size_t j = 0; j < gm.n_rows; j++)
  //   {
  //     internal_data.push_back(gm(j, i));
  //   }
  // }

  // std::unique_ptr<ranger::Data> internal_data_ptr { };
  // internal_data_ptr = ranger::make_unique<SeidrForestData>();
  // internal_data_ptr->load_seidr(dependent_variable_name, gm_x, gm_y, var);

  init(M, load(dependent_variable_name, gm_x, gm_y, var), mtry,
       //output_prefix, num_tree, seed, threads, importance mode
       "ranger_out_s", ntree, 314159265, 1, I,
       //min node size (0=auto), status variable, prediction mode,
       // sample with replacement
       min_node_size, false, true,
       //categorical vars, save memory, split rule (1 = default/variance)
       catvars, false, S,
       //predict all, sample fraction, alpha, min prop, holdout
       false, sample_fraction_vector, alpha, minprop, false,
       //prediciton mode (1=response), number of random splits
       P, 1, false, 0);
  // Catch random log messages in a string sink
  verbose_out = &ostream;
}



void genie3(const arma::mat& gm,
            const std::vector<std::string>& genes,
            const std::vector<arma::uword>& pred,
            const std::string& tmpdir,
            const uint64_t& ntree,
            const uint64_t& mtry,
            const uint64_t& min_node_size,
            const double& alpha,
            const double& minprop)
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  std::string tmpfile = tempfile(tmpdir);
  std::ofstream ofs(tmpfile.c_str(), std::ios::out);
  if (! ofs)
  {
    throw std::runtime_error("Could not open temp file: " + tmpfile);
  }
  //SeidrForestData data(gm, genes);
  #pragma omp parallel for
  for (arma::uword i = 0; i < pred.size(); i++)
  {

    #pragma omp critical
    {
      log << "Gene: " << pred[i] << '\n';
      log.send(LOG_INFO);
    }

    std::vector<std::string> dep_var;
    std::vector<std::string> var;

    arma::uvec ranger_ix(gm.n_cols - 1);
    arma::uvec ranger_iy(1);

    ranger_iy(0) = pred[i];
    dep_var.push_back(genes[pred[i]]);

    arma::uword j_ = 0;
    for (arma::uword j = 0; j < gm.n_cols; j++)
    {
      if (j != pred[i])
      {
        ranger_ix(j_) = j;
        var.push_back(genes[j]);
        j_++;
      }
    }
    const arma::mat& ranger_mx = gm.cols(ranger_ix);
    const arma::mat& ranger_my = gm.cols(ranger_iy);

    std::vector<std::string> catvars;
    std::ostream nullstream(0);
    SeidrForest forest(dep_var, ranger_mx, ranger_my, var, ntree, mtry,
                       min_node_size, alpha, minprop, nullstream);

    // forest->init(genes[pred[i]], MEM_DOUBLE, data, mtry, "ranger_out_s",
    //              ntree, 314159265, 1, IMP_GINI, min_node_size, "", false,
    //              true, catvars, false, LOGRANK, false, 1, alpha, minprop,
    //              false, RESPONSE, 1);

    forest.run(false, false);
    std::vector<double> varimp = forest.getVariableImportance();
    #pragma omp critical
    {
      ofs << pred[i] << '\n';

      arma::uword j = pred[i];

      if (j == 0)
      {
        ofs << 0 << '\t';
        for (arma::uword k = 0; k < varimp.size(); k++)
        {
          ofs << varimp[k]
              << (k == varimp.size() - 1 ? '\n' : '\t');
        }
      }
      else if (j == varimp.size())
      {
        for (arma::uword k = 0; k < varimp.size(); k++)
        {
          ofs << varimp[k] << '\t';
        }
        ofs << 0 << '\n';
      }
      else
      {
        for (arma::uword k = 0; k < varimp.size(); k++)
        {
          if (j == k)
            ofs << 0 << '\t';
          ofs << varimp[k]
              << (k == varimp.size() - 1 ? '\n' : '\t');
        }
      }
    }
  }
  ofs.close();
}

void genie3_full(const arma::mat& gm,
                 const std::vector<std::string>& genes,
                 const seidr_genie3_param_t& param)
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());
  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> uvec;
  for (uint64_t i = 0; i < gm.n_cols; i++)
  {
    uvec.push_back(i);
  }

  seidr_mpi_genie3 mpi(param.bs, gm, uvec, genes, param.tempdir,
                       param.outfile);

  mpi.set_ntree(param.ntree);
  mpi.set_min_node_size(param.min_node_size);
  mpi.set_alpha(param.alpha);
  mpi.set_minprop(param.minprop);
  mpi.set_mtry(param.mtry);

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

void genie3_partial(const arma::mat& gm,
                    const std::vector<std::string>& genes,
                    const std::vector<std::string>& targets,
                    const seidr_genie3_param_t& param)
{
  seidr_mpi_logger log(LOG_NAME"@" + mpi_get_host());

  fs::path p_out(param.outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<uint64_t> positions;
  for (uint64_t i = 0; i < targets.size(); i++)
  {
    uint64_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if (pos >= genes.size())
    {
      log << "Gene " << targets[i] << " was not found in the expression set "
          << "and will therefore not be considered. "
          << "Please check that your expression set and "
          << "its column names (gene file) contain an entry for "
          << targets[i] << ".\n";
      log.log(LOG_WARN);
    }
    else
    {
      positions.push_back(pos);
    }
  }

  seidr_mpi_genie3 mpi(param.bs, gm, positions, genes, param.tempdir,
                       param.outfile);
  mpi.set_ntree(param.ntree);
  mpi.set_min_node_size(param.min_node_size);
  mpi.set_alpha(param.alpha);
  mpi.set_minprop(param.minprop);
  mpi.set_mtry(param.mtry);
  mpi.set_targetd(true);

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
