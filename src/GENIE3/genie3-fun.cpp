// Seidr
#include <common.h>
#include <genie3-fun.h>
#include <mpims.h>
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
#include <DataDouble.h>
#include <utility.h>
#include <ForestRegression.h>
#include <globals.h>
#include <map>
#include <sstream>

namespace fs = boost::filesystem;
typedef std::pair<std::streampos, std::string> file_index_t;
typedef std::map<seidr_uword_t, file_index_t> result_file_map_t;

class seidr_mpi_genie3 : public seidr_mpi {
public:
  void entrypoint();
  void finalize();
  seidr_mpi_genie3() : seidr_mpi() {}
  seidr_mpi_genie3(unsigned long bs) : seidr_mpi(bs) {}
  void set_ntree(seidr_uword_t x){_ntree = x;}
  void set_mtry(seidr_uword_t x){_mtry = x;}
  void set_min_node_size(seidr_uword_t x){_min_node_size = x;}
  void set_alpha(double x){_alpha = x;}
  void set_minprop(double x){_minprop = x;}
private:
  seidr_uword_t _ntree;
  seidr_uword_t _mtry;
  seidr_uword_t _min_node_size;
  double _alpha;
  double _minprop;
};

void seidr_mpi_genie3::entrypoint()
{
  seidr_mpi_logger log;
  std::string tmpfile = _outfilebase + "/.seidr_tmp_genie3/MPIthread_" +
    std::to_string(_id) + "_" + std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);
  
  std::vector<arma::uword> uvec;
  for(auto i : _indices)
    uvec.push_back(i);
  genie3(_data, _genes, uvec, tmpfile, _id, _ntree, _mtry,
         _min_node_size, _alpha, _minprop);
  announce_ready();
}

void seidr_mpi_genie3::finalize()
{
  seidr_mpi_logger log;
  result_file_map_t rmap;
  if(_id == 0)
    {
      seidr_mpi_logger log;
      log << "Merging tmp files and cleaning up.\n";
      log.send(LOG_INFO);
      std::vector<fs::path> files;
      std::string tmpdir = _outfilebase + "/.seidr_tmp_genie3";
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

class SeidrForestData : public DataDouble {
public:
  SeidrForestData(arma::mat& gm, std::vector<std::string>& genes);
};

SeidrForestData::SeidrForestData(arma::mat& gm,
                                 std::vector<std::string>& genes)
{
  num_rows = gm.n_rows;
  externalData = false;
  variable_names = genes;
  num_cols = variable_names.size();
  num_cols_no_snp = num_cols;
  reserveMemory();
  bool error = false;
  for(arma::uword i = 0; i < gm.n_rows; i++)
    for(arma::uword j = 0; j < gm.n_cols; j++)
      set(j, i, gm(i, j), error);
}

class SeidrForest : public ForestRegression {
public:
  SeidrForest(std::string dependent_variable_name, Data* input_data,
              size_t ntree, size_t mtry,
              size_t min_node_size, double alpha, double minprop,
              std::ostream& ostream);
};

SeidrForest::SeidrForest(std::string dependent_variable_name,
                         Data* input_data, size_t ntree, size_t mtry,
                         size_t min_node_size, double alpha, double minprop,
                         std::ostream& ostream){
  std::vector<std::string> catvars;
  MemoryMode M = MEM_DOUBLE;
  ImportanceMode I = IMP_GINI;
  SplitRule S = LOGRANK;
  PredictionType P = RESPONSE;
  
  init(dependent_variable_name, M, input_data, mtry,
       //output_prefix, num_tree, seed, threads, importance mode
       "ranger_out_s", ntree, 314159265, 1, I,
       //min node size (0=auto), status variable, prediction mode,
       // sample with replacement
       min_node_size, "", false, true,
       //categorical vars, save memory, split rule (1 = default/variance)
       catvars, false, S,
       //predict all, sample fraction, alpha, min prop, holdout
       false, 1, alpha, minprop, false,
       //prediciton mode (1=response), number of random splits
       P, 1);
  // Catch random log messages in a string sink
  verbose_out = &ostream;
}



void genie3(arma::mat gm, std::vector<std::string> genes,
            std::vector<arma::uword> pred, std::string outfile,
            int thread_id,
            seidr_uword_t ntree, seidr_uword_t mtry,
            size_t min_node_size, double alpha, double minprop)
{
  seidr_mpi_logger log;
  std::ofstream ofs(outfile.c_str(), std::ios::out);
  if(! ofs)
    throw std::runtime_error("Could not open temp file: " + outfile);
  //SeidrForestData data(gm, genes);

  for(arma::uword i = 0; i < pred.size(); i++)
    {

      log << "Gene: " << pred[i] << '\n';
      log.send(LOG_INFO);
      Data* data = new DataDouble(gm.memptr(), genes, gm.n_rows, gm.n_cols);
      //Forest* forest = new ForestRegression;
      //forest->verbose_out=&std::cout;
      std::vector<std::string> catvars;
      std::ostream nullstream(0);
      SeidrForest forest(genes[pred[i]], data, ntree, mtry, min_node_size,
                         alpha, minprop, nullstream);

      // forest->init(genes[pred[i]], MEM_DOUBLE, data, mtry, "ranger_out_s", 
      //              ntree, 314159265, 1, IMP_GINI, min_node_size, "", false, 
      //              true, catvars, false, LOGRANK, false, 1, alpha, minprop,
      //              false, RESPONSE, 1);

      forest.run(false);
      std::vector<double> varimp = forest.getVariableImportance();
      ofs << pred[i] << '\n';

      arma::uword j = pred[i];
      
      if(j == 0)
        {
          ofs << 0 << '\t';
          for(arma::uword k = 0; k < varimp.size();k++)
            {
              ofs << varimp[k]
                  << (k == varimp.size() - 1 ? '\n' : '\t');
            }
        }
      else if(j == varimp.size())
        {
          for(arma::uword k = 0; k < varimp.size();k++)
            {
              ofs << varimp[k] << '\t';
            }
          ofs << 0 << '\n';
        }
      else
        {
          for(arma::uword k = 0; k < varimp.size(); k++)
            {
              if(j == k)
                ofs << 0 << '\t'; 
              ofs << varimp[k]
                  << (k == varimp.size() - 1 ? '\n' : '\t');
            }
        }
      //delete forest; 
      delete data;
    }
    
}
  
void genie3_full(arma::mat gm, std::vector<std::string> genes,
                 arma::uword bs, std::string outfile,
                 seidr_uword_t ntree, seidr_uword_t mtry,
                 size_t min_node_size, double alpha, double minprop)
{
  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<unsigned long> uvec;
  for(unsigned long i = 0; i < gm.n_cols; i++)
    uvec.push_back(i);
  seidr_mpi_genie3 mpi(bs);

  mpi.set_data(gm);
  mpi.set_genes(genes);
  mpi.set_indices(uvec);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_ntree(ntree);
  mpi.set_min_node_size(min_node_size);
  mpi.set_alpha(alpha);
  mpi.set_minprop(minprop);
  mpi.set_mtry(mtry);

  mpi.exec();
}
  
void genie3_partial(arma::mat gm, std::vector<std::string> genes, size_t bs,
                    std::vector<std::string> targets, std::string outfile,
                    seidr_uword_t ntree, seidr_uword_t mtry,
                    size_t min_node_size, double alpha, double minprop)
{
  seidr_mpi_logger log;

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());

  std::vector<size_t> positions;
  for(size_t i = 0; i < targets.size(); i++){
    size_t pos = find(genes.begin(), genes.end(), targets[i]) - genes.begin();
    if(pos >= genes.size()){
      log << "Gene " << targets[i] << " was not found in the expression set " <<
        "and will therefore not be considered. Please check that your expression set and " <<
        "its column names (gene file) contain an entry for " << targets[i] << ".\n";
      log.log(LOG_WARN);
    } else {
      positions.push_back(pos);
    }
  }

  seidr_mpi_genie3 mpi(bs);
  mpi.set_data(gm);
  mpi.set_genes(genes);
  mpi.set_indices(positions);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_ntree(ntree);
  mpi.set_mtry(mtry);
  mpi.set_min_node_size(min_node_size);
  mpi.set_alpha(alpha);
  mpi.set_minprop(minprop);

  mpi.exec();
}
