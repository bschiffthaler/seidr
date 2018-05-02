// Seidr
#include <common.h>
#include <tiglm.h>
#include <glmnet2.h>
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

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.2,1.0);

namespace fs = boost::filesystem;

class seidr_mpi_tigress : public seidr_mpi {
public:
  void entrypoint();
  void finalize();
  seidr_mpi_tigress() : seidr_mpi() {}
  seidr_mpi_tigress(unsigned long bs) : seidr_mpi(bs) {}
  void set_nboot(seidr_uword_t x){_nboot = x;}
  void set_fmin(double x){_fmin = x;}
  void set_nsteps(seidr_uword_t x){_nsteps = x;}
  seidr_uword_t get_nboot(){return _nboot;}
  double get_fmin(){return _fmin;}
  seidr_uword_t set_nsteps(){return _nsteps;}
private:
  seidr_uword_t _nboot;
  seidr_uword_t _nsteps;
  double _fmin;
};

void seidr_mpi_tigress::entrypoint()
{
  seidr_mpi_logger log;
  std::string tmpfile = _outfilebase + "/.seidr_tmp_tigress/MPIthread_" +
    std::to_string(_id) + "_" + std::to_string(_indices[0]) + ".txt";

  log << "Using tempfile '" << tmpfile << "'\n";
  log.send(LOG_INFO);
  
  std::vector<arma::uword> uvec;
  for(auto i : _indices)
    uvec.push_back(i);
  tiglm(_data, _genes, uvec, tmpfile, _id, _nboot, _fmin, _nsteps);
  announce_ready();
}

void seidr_mpi_tigress::finalize()
{
  if(_id == 0)
    {
      seidr_mpi_logger log;
      log << "Merging tmp files and cleaning up.\n";
      log.send(LOG_INFO);
      std::vector<fs::path> files;
      fs::path p_tmp(_outfilebase + "/.seidr_tmp_tigress");
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
            ofs << l << '\n';
          fs::remove(p);
        }
      fs::remove(p_tmp);
      log << "Finished.\n";
      log.send(LOG_INFO);
    }
  double now = MPI_Wtime();
  while(MPI_Wtime() - now < 10)
    {
      check_logs();
    }
}

/* 
In order to do stability selection we need to split the samples
in half at random. This function takes a sample size and outputs
a std::vector containing two arma::uvec vectors which hold the
randomized indices of the matrix rows.
*/
std::vector<arma::uvec> shuffle(unsigned int n){
  arma::uvec samples(n);
  unsigned int nA;
  unsigned int nB;

  // Determine if sample size is even or uneven to get the right
  // size for the uvec objects
  if( n % 2 == 0 ) {
    nA = n/2;
    nB = n/2;
  } else {
    unsigned int nT = n+1;
    nA = nT/2;
    nB = nT/2 - 1;
  }
  arma::uvec resA(nA);
  arma::uvec resB(nB);
  std::vector<arma::uvec> final;

  for (unsigned int i = 0; i < n; i++) samples[i]=i;

  arma::uvec samplesS = arma::shuffle(samples);

  for(unsigned int i = 0; i < n; i++){
    if(i < nA) resA[i] = samplesS[i];
    if(i >= nA) resB[i-nA] = samplesS[i];
  }

  final.push_back(resA);
  final.push_back(resB);
  return final;
  
}


void tiglm(arma::mat geneMatrix, std::vector<std::string> genes,
	   std::vector<arma::uword> pred, std::string outfile,
           int thread_id, seidr_uword_t boot, double fmin,
           seidr_uword_t nsteps){

  seidr_mpi_logger log;
  size_t ng = geneMatrix.n_cols - 1;
  size_t ns = geneMatrix.n_rows;
  double ns_d = geneMatrix.n_rows;
  
  size_t n_good = ceil(sqrt(ns_d / 2));

  log << "Minimal required number of useful samples: "
      << n_good << '\n';;
  log.send(LOG_DEBUG);
  std::ofstream ofs(outfile.c_str(), std::ios::out);

  if(! ofs)
    throw std::runtime_error("Could not open temp file: " + outfile);
  
  for(size_t i = 0; i < pred.size(); i++){
    log << "Gene: " << pred[i] << '\n';
    log.send(LOG_DEBUG);
    clock_t begin = clock();
    arma::uvec indices = get_i(pred[i],ng + 1);
    arma::mat cum = arma::zeros<arma::mat>(ng, nsteps);

    arma::mat response = geneMatrix.cols(indices);
    arma::vec predictor = geneMatrix.col(pred[i]);

    size_t iter  = 0;
    size_t nboot = 0;
    /* This int is incremented if the random shuffle produced a vector
       with all 0s. If that happens too often, we skip the gene and
       return no interations, as the other case could result in an
       endless loop or a non-converging lasso. */

    arma::vec randV(ng,1);

    std::seed_seq sseq{3, 1, 4, 1, 5, 9, 2, 6, 5};
    generator.seed(sseq);
    arma::arma_rng::set_seed(314159265);
    
    while(iter < boot){

      // Shuffle the samples
      std::vector<arma::uvec > split = shuffle(ns);
      arma::vec yA = predictor(split[0]);
      arma::vec yB = predictor(split[1]);

      arma::uvec syA = arma::find_unique(yA);
      arma::uvec syB = arma::find_unique(yB);
      
      // If any of the vectors have less than two usable values, try again
      if(syA.n_elem < n_good || syB.n_elem < n_good) {
	nboot++;
	if(nboot > 1000)
	  {
	    log << genes[pred[i]] << " encountered too many reshuffles... skipping.\n";
            log.send(LOG_WARN);
	    cum.fill(0);
	    iter = boot;
	  }
        continue;
      }

      randV.imbue( [&]() { return distribution(generator); } );
      arma::mat rw = response.t();
      for(unsigned int y = 0; y < ns; y++){
        rw.col(y) % randV;
      }
      rw = rw.t();

      arma::mat xA = rw.rows(split[0]);
      arma::mat xB = rw.rows(split[1]);

      glm pathA = glmnet(xA, yA, nsteps, fmin);
      glm pathB = glmnet(xB, yB, nsteps, fmin);

      if(pathA.success && pathB.success)
        {
          cum += arma::sign(arma::abs(pathA.beta));
          cum += arma::sign(arma::abs(pathB.beta));
          iter++;
        }
    }
    double d_nsteps = 0;
    try
      {
        d_nsteps = boost::numeric_cast<double>(nsteps);
      }
    catch(boost::numeric::bad_numeric_cast& e)
      {
        throw std::runtime_error(e.what());
      }
    
    cum = arma::cumsum(cum.t())/(boot*2)/d_nsteps;

    for(arma::uword s = 0; s < indices.n_elem; s++){
      arma::rowvec fin = cum.row(nsteps - 1);
      ofs << genes[pred[i]] << "\t"
          << genes[indices[s]] << "\t"
          << fin(s) << std::endl;
    }

    clock_t end = clock();
    log << "Processed gene " << pred[i] << " in "
        << double(end - begin) / CLOCKS_PER_SEC << " s"
        << '\n';
    log.send(LOG_INFO);
  }
}

void tiglm_full(arma::mat GM, std::vector<std::string> genes, size_t bs,
		std::string outfile, seidr_uword_t boot, double fmin,
                seidr_uword_t nsteps){

  fs::path p_out(outfile);
  p_out = fs::absolute(p_out);

  fs::path d_out(p_out.parent_path());
  
  std::vector<unsigned long> uvec;
  for(unsigned long i = 0; i < GM.n_cols; i++)
    uvec.push_back(i);
  seidr_mpi_tigress mpi(bs);
  
  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(uvec);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_nboot(boot);
  mpi.set_fmin(fmin);
  mpi.set_nsteps(nsteps);
  
  mpi.exec();
} 

void tiglm_partial(arma::mat GM, std::vector<std::string> genes, size_t bs,
		   std::vector<std::string> targets, std::string outfile,
                   seidr_uword_t boot, double fmin,
                   seidr_uword_t nsteps){

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

  seidr_mpi_tigress mpi(bs);
  mpi.set_data(GM);
  mpi.set_genes(genes);
  mpi.set_indices(positions);
  mpi.set_outfilebase(d_out.string());
  mpi.set_outfile(p_out.string());
  mpi.set_nboot(boot);
  mpi.set_fmin(fmin);
  mpi.set_nsteps(nsteps);
  
  mpi.exec();
}
