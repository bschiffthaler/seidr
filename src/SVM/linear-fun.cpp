// Seidr
#include <common.h>
#include <linear-fun.h>
#include <mpims.h>
// External
#include <iostream>
#include <random>
#include <string>
#include <linear.h>
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

void print_null(const char *s) {}

class seidr_mpi_svm : public seidr_mpi {
public:
  void entrypoint();
  void finalize();
  seidr_mpi_svm() : seidr_mpi() {}
  seidr_mpi_svm(unsigned long bs) : seidr_mpi(bs) {}
  void set_param(parameter& p) {_param = p;}
  void set_min_sample_size(arma::uword x) {_min_sample_size = x;}
  void set_max_sample_size(arma::uword x) {_max_sample_size = x;}
  void set_predictor_sample_size_min(arma::uword x) {_predictor_sample_size_min = x;}
  void set_predictor_sample_size_max(arma::uword x) {_predictor_sample_size_max = x;}
  void set_ensemble_size(arma::uword x) {_ensemble_size = x;}
private:
  parameter _param;
  arma::uword _min_sample_size;
  arma::uword _max_sample_size;
  arma::uword _predictor_sample_size_min;
  arma::uword _predictor_sample_size_max;
  arma::uword _ensemble_size;
};

void seidr_mpi_svm::entrypoint()
{
  seidr_mpi_logger log;
  std::string tmpfile = _outfilebase + "/.seidr_tmp_llr/MPIthread_" +
    std::to_string(_id) + "_" + std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);
  
  std::vector<arma::uword> uvec;
  for(auto i : _indices)
    uvec.push_back(i);
  svm(_data, _genes, uvec, tmpfile, _id, _param, _min_sample_size,
      _max_sample_size, _predictor_sample_size_min,
      _predictor_sample_size_max, _ensemble_size);
  announce_ready();
}

void seidr_mpi_svm::finalize()
{
  result_file_map_t rmap;
  if(_id == 0)
    {
      seidr_mpi_logger log;
      log << "Merging tmp files and cleaning up.\n";
      log.send(LOG_INFO);
      std::vector<fs::path> files;
      std::string tmpdir = _outfilebase + "/.seidr_tmp_llr";
      fs::path p_tmp(tmpdir);
      for(auto it = fs::directory_iterator(p_tmp);
          it != fs::directory_iterator(); it++)
        {
          if( fs::is_regular_file( it->path() ) )
            files.push_back( (*it).path() );
        }
      std::ofstream ofs(_outfile);

      for(fs::path& p : files)
        {
          std::ifstream ifs(p.string().c_str());
          std::string l;
          while(std::getline(ifs,l))
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
      for(; it != rmap.end(); it++)
        {
          ifs.seekg(gene_index);
          std::string l;
          std::getline(ifs, l);
          ofs << l << '\n';

          auto nx = std::next(it);
          if(nx->second.second != file_path)
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
      while(MPI_Wtime() - now < 10)
        {
          check_logs();
        }
    }
}

void svm(arma::mat geneMatrix, std::vector<std::string> genes,
         std::vector<arma::uword> uvec, std::string outfile,
         int thread_id, parameter param,
         arma::uword min_sample_size, arma::uword max_sample_size,
         arma::uword predictor_sample_size_min,
         arma::uword predictor_sample_size_max,
         arma::uword ensemble_size){
  
  seidr_mpi_logger log;  

  std::ofstream ofs(outfile.c_str(), std::ios::out);
  set_print_string_function(print_null);
  if(! ofs)
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
  catch(boost::numeric::bad_numeric_cast& e)
    {
      throw std::runtime_error(e.what());
    }
  

  arma::uvec samples(geneMatrix.n_rows);
  for(arma::uword i = 0; i < geneMatrix.n_rows; i++)
    {
      samples(i) = i;
    }

  for(auto& target : uvec)
    {
      log << "Started gene: " << genes[target] << ".\n";
      log.send(LOG_INFO);
      arma::uvec pred(geneMatrix.n_cols - 1);
      arma::uword j = 0;
      for(arma::uword i = 0; i < geneMatrix.n_cols; i++)
        {
          if(i != target)
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
      for(arma::uword boot = 0; boot < ensemble_size; boot++)
        {
          // Randomly select a number of predictor genes in range
          // predictor_sample_size_min to predictor_sample_size_max
          pred = arma::shuffle(pred);
          arma::uword pred_size = predictor_gen(gen);
          arma::uvec pred_sub(pred_size);
          for(arma::uword i = 0; i < pred_size; i++)
            {
              pred_sub(i) = pred(i);
            }

          // Randomly select a number of samples in range
          // min_sample_size to max_sample_size
          samples = arma::shuffle(samples);
          arma::uword sample_size = sample_gen(gen);
          arma::uvec sample_sub(sample_size);
          for(arma::uword i = 0; i < sample_size; i++)
            {
              sample_sub(i) = samples(i);
            }

          // Create subset matrix
          arma::mat pred_mat = geneMatrix.submat(sample_sub, pred_sub);

          // Set up SVM problem struct sparse node notation
          feature_node ** nodes = SEIDR_MALLOC(feature_node*, pred_mat.n_rows);

          for(arma::uword i = 0; i < pred_mat.n_rows; i++)
            {
              nodes[i] = SEIDR_MALLOC(feature_node, pred_mat.n_cols + 2);
            }
  
          for(arma::uword i = 0; i < pred_mat.n_rows; i++)
            {
              arma::uword t = 0;
              for(arma::uword j = 0; j < pred_mat.n_cols; j++)
                {
                  int x = 0;
                  // Safely cast uword to int
                  try
                    {
                      x = boost::numeric_cast<int>(j + 1);
                    }
                  catch(boost::numeric::bad_numeric_cast& e)
                    {
                      throw std::runtime_error(e.what());
                    }
                  double y = pred_mat(i,j);
                  // Do not add zeros to the problem
                  if(! almost_equal(y, 0))
                    {
                      nodes[i][t] = feature_node{x, y};
                      t++;
                    }
                }
              int bias_index = 0;
              try
                {
                  bias_index = boost::numeric_cast<int>(pred_mat.n_cols + 1);
                }
              catch(boost::numeric::bad_numeric_cast& e)
                {
                  throw std::runtime_error(e.what());
                }
              nodes[i][t] = feature_node{bias_index, 1};
              t++;
              // The end of the nodes is signified by a node with index -1
              nodes[i][t] = feature_node{-1, 0};
              //std::cout << "Final t:" << t << '\n';
            }

          // Get target expression and make sure we have the same samples
          // as in the predictor expression
          arma::vec target_exp = geneMatrix.col(target);
          target_exp = target_exp.elem(sample_sub);

          int nrows = 0;
          int nfeats = 0;
          try
            {
              nrows = boost::numeric_cast<int>(pred_mat.n_rows);
              // We add 1 to the number of features for the bias value
              nfeats = boost::numeric_cast<int>(pred_mat.n_cols) + 1;
            }
          catch(boost::numeric::bad_numeric_cast& e)
            {
              throw std::runtime_error(e.what());
            }
      
          problem prob{nrows, nfeats, target_exp.memptr(), nodes, 1.0};

          // Check parameters of the model
          const char * error_msg = check_parameter(&prob, &param);
          if(error_msg)
            {
              std::string msg(error_msg);
              throw std::runtime_error(msg);
            }

          model * mod = train(&prob, &param);

          // Extract support vector coefficients and indices
          arma::mat coefs(1, nfeats - 1);
          for(int i = 0; i < nfeats - 1; i++)
            {
              coefs(0, i) = mod->w[i];
            }

          //coefs.submat(0,0,0,5).print();
          
          /*arma::uvec sv_ind(mod->l);
          int nr_sv = svm_get_nr_sv(model);
          int * sv_indices = new int[nr_sv];
          svm_get_sv_indices(model, sv_indices);
          for(int i = 0; i < nr_sv; i++)
            {
              sv_ind[i] = sv_indices[i] - 1;
              }*/

          // Get support vectors from the predictor matrix
          //arma::mat svs = pred_mat.rows(sv_ind);

          //arma::mat pred_transp = pred_mat.t();
          // Calculate absolute value of weights and rank
          arma::mat weights = coefs;
          weights = arma::abs(weights);
          arma::uvec sort_index = arma::sort_index(weights, "descend");

          //weights.tail_cols(5).print();
          
          // Check which weights were zero
          std::vector<arma::uword> zeros;
          for(arma::uword i = 0; i < weights.n_elem; i++)
            {
              if(almost_equal(weights(i), 0))
                {
                  zeros.push_back(i);
                }
            }
          arma::uvec zero_ind = arma::conv_to<arma::uvec>::from(zeros);
          
          // Zero out weights vector and set top 20th values to 1
          weights.zeros();
          arma::uword rt = pred_mat.n_cols / 20;
          for(arma::uword i = 0; i <= rt; i++)
            {
              weights(sort_index(i)) = 1;
            }
          // Reset those weights which were zero before ranking
          weights.elem(zero_ind).zeros();
          // Update return vector with new model data
          ret.elem(pred_sub) += weights;

          // Free memory of SVM problem and model
          free_and_destroy_model(&mod);
          //delete[] sv_indices;
          //std::cout << "Freeing " << pred_mat.n_rows << " nodes*\n";
          for(arma::uword i = 0; i < pred_mat.n_rows; i++)
            {
              //std::cout << i;
              free(nodes[i]);
              //std::cout << "..";
            }
          //std::cout << "\nFreed nodes*\n";
          free(nodes);
          //std::cout << "Freed nodes**\n";
        }

      ofs << target << '\n';
      for(seidr_uword_t i = 0; i < ret.size(); i++)
        {
          ofs << ret[i] / fensemble_size
              << (i == ret.size() - 1 ? '\n' : '\t');
        }
    }
  destroy_param(&param);
  
}

void svm_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
              std::string outfile, parameter param,
              arma::uword min_sample_size, arma::uword max_sample_size,
              arma::uword predictor_sample_size_min,
              arma::uword predictor_sample_size_max,
              arma::uword ensemble_size){

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());
  
  std::vector<unsigned long> uvec;
  for(unsigned long i = 0; i < GM.n_cols; i++)
    uvec.push_back(i);
  seidr_mpi_svm mpi(bs);
  
  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(uvec);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_param(param);
  mpi.set_min_sample_size(min_sample_size);
  mpi.set_max_sample_size(max_sample_size);
  mpi.set_predictor_sample_size_min(predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(predictor_sample_size_max);
  mpi.set_ensemble_size(ensemble_size);
  
  mpi.exec();
} 

void svm_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
                 std::vector<std::string> targets, std::string outfile,
                 parameter param, arma::uword min_sample_size,
                 arma::uword max_sample_size,
                 arma::uword predictor_sample_size_min,
                 arma::uword predictor_sample_size_max,
                 arma::uword ensemble_size){

  seidr_mpi_logger log;
  
  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<size_t> positions;
  for(size_t i = 0; i < targets.size(); i++){
    size_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if(pos >= genes.size()){
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
  
  seidr_mpi_svm mpi(bs);
  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(positions);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_param(param);
  mpi.set_min_sample_size(min_sample_size);
  mpi.set_max_sample_size(max_sample_size);
  mpi.set_predictor_sample_size_min(predictor_sample_size_min);
  mpi.set_predictor_sample_size_max(predictor_sample_size_max);
  mpi.set_ensemble_size(ensemble_size);
  mpi.exec();  
}
