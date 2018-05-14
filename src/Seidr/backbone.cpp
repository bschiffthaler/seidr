// Seir
#include <common.h>
#include <RankFile.h>
#include <gzstream.h>
#include <fs.h>
#include <BSlogger.h>
#include <lower_triangular.h>
// External
#include <vector>
#include <string>
#include <fstream>
#include <tclap/CmdLine.h>
#include <stdexcept>
#include <cerrno>
#include <map>
#include <boost/lexical_cast.hpp>
#include <armadillo>
#include <gsl/gsl_integration.h>

using boost::lexical_cast;

// Integration function
double f (double x, void * params) {
  uint64_t k = *(uint64_t *) params;
  double r = pow(1-x,k-2);
  return r;
}

int backbone(int argc, char * argv[])
{

  logger log(std::cerr, "backbone");
  
  // Variables used by the function
  std::string infs;
  seidr_score_t alpha;
  std::string out_file;
  bool reverse;
  
  // We ignore the first argument, the function name
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " backbone";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
    {
      // Add arguments from the command line
      TCLAP::CmdLine cmd("Filter the network by its multiscale backbone (Serrano et al. 2009)", ' ', version);
      TCLAP::ValueArg<std::string> 
      arg_outfile("o", "outfile", "Where to write the output file.", false,
                                              "backbone.bin", "string");
      cmd.add(arg_outfile);
      TCLAP::ValueArg<std::string> arg_input_file("i","infile","Input file",true,"","string");
      cmd.add(arg_input_file);
      TCLAP::ValueArg<seidr_score_t> arg_alpha("a","alpha","Alpha cutoff for filtering",false,0.05,"float");
      cmd.add(arg_alpha);
      TCLAP::SwitchArg switch_reverse("r","reverse","Reverse scores (lower is better)",cmd,false);
      
      // Parse arguments
      cmd.parse(argc, args);
      infs = arg_input_file.getValue();
      out_file = arg_outfile.getValue();
      alpha = arg_alpha.getValue();
      reverse = switch_reverse.getValue();
    }
  catch(TCLAP::ArgException& except)
    {
      log(LOG_ERR) << except.what() << '\n';
      return EINVAL;
    }

  try // FS checks
    {
    }
  catch (std::runtime_error& except)
    {
      log << except.what() << '\n';
      return errno;
    }

  std::vector<std::string> reference_genes;
  std::vector<std::string> algorithms;
  
  std::string cif = infs.c_str();
  igzstream rf(cif.c_str());
  
  cereal::BinaryInputArchive iarchive(rf);

  LowerTriangular< RankFileEdge<seidr_score_t> > lt_mat;
  
  log(LOG_INFO) << "Reading " << infs << '\n';
  size_t current_i = 0;
  size_t current_j = 0;
  size_t current_n = 0;
  seidr_score_t max = std::numeric_limits<seidr_score_t>::infinity();
  seidr_score_t min = -max;
  try
    {
      RankFileHeader h;
      iarchive(h);
	
      reference_genes = h.get_node_names();
      algorithms = h.get_method_names();

      lt_mat.reset(h.get_nn());
      for(uint64_t i = 1; i < lt_mat.size(); i++)
        {
          for(uint64_t j = 0; j < i; j++)
            {
              seidr_score_t n = std::numeric_limits<seidr_score_t>::quiet_NaN();
              RankFileEdge<seidr_score_t>& e = lt_mat.at(i, j);
              e.set_r(n);
            }
        }

      
      bool first = true;
      
      for(size_t i = 0; i < h.get_ne(); i++)
        {
          current_n = i;
          RankFileEdge<seidr_score_t> e;
          iarchive(e);
          size_t ei = e.get_i();
          size_t ej = e.get_j();
          seidr_score_t es = e.get_r();
          current_i = ei;
          current_j = ej;
          // Lower triangular order
          // This should be certain coming from el2bin. I'm just paranoid
          if(ei < ej) 
            {
              size_t _ei = ej;
              ej = ei;
              ei = _ei;
              if(e.get_d() == 2)
                {
                  e.set_d(1);
                }
              else if(e.get_d() == 1)
                {
                  e.set_d(2);
                }
            }
          if(first)
            {
              min = es;
              max = es;
              first = false;
            }
          else
            {
              if(es < min)
                min = es;
              if(es > max)
                max = es;
            }
          lt_mat.at(ei, ej) = e;
        }
    }
  catch(std::exception& except)
    {
      log(LOG_ERR) << except.what() << " at index "
                   << current_i << ',' << current_j
                   << " after reading " << current_n
                   <<" items\n";
      return errno;
    }


  log(LOG_INFO) << "Done reading, have values in range ["
                << min << ',' << max << "]\n";
  
  // Substract scores from maximum to switch order in
  // reverse mode
  if(reverse)
    {
      for(uint64_t i = 1; i < lt_mat.size(); i++)
        {
          for(uint64_t j = 0; j < i; j++)
            {
              auto& v = lt_mat.at(i, j);
              if(! std::isnan( v.get_r() ) )
                {
                  v.set_r( max - v.get_r() );
                }
            }
        }
    }
  
  for(uint64_t i = 0; i < lt_mat.size(); i++)
    {
      seidr_score_t ws = 0; //sum of weights
      uint64_t k = 0; //number of edges

      log(LOG_INFO) << "Testing node " << reference_genes[i] << '\n';
      
      // First pass to get parameters
      for(uint64_t j = 0; j < lt_mat.size(); j++)
        {
          if(i == j) continue;
          // Override indices i and j internally here to ensure
          // lower triangular sorting
          uint64_t _i,_j;
          if(i < j)
            {
              _i = j;
              _j = i;
            }
          else
            {
              _i = i;
              _j = j;
            }
          auto& v = lt_mat.at(_i, _j);
          if(! std::isnan( v.get_r() ) )
            {
              k++;
              ws += v.get_r();
            }
        }

      
      // Second pass to get p-value
      if(k > 1)
        {
          for(uint64_t j = 0; j < lt_mat.size(); j++)
            {
              if(i == j) continue;
              uint64_t _i,_j;
              if(i < j)
                {
                  _i = j;
                  _j = i;
                }
              else
                {
                  _i = i;
                  _j = j;
                }
              auto& v = lt_mat.at(_i, _j);
              if(std::isnan( v.get_r() ) )
                continue;

              auto w = v.get_r(); // edge weight

              double pij = w/ws;

              double res, err;
              gsl_integration_workspace * work
                = gsl_integration_workspace_alloc (k * 2);
              gsl_function F;
              F.function = &f;
              F.params = &k;

              gsl_integration_qags (&F, 0, pij, 0, 1e-7, k * 2 - 1,
                                    work, &res, &err);

              double aij = 1 - (k-1) * res;
              v.push_back_w(aij);
              gsl_integration_workspace_free (work);
            }
        }
    }


  
  try
    {
     ogzstream rf;
     rf.open(out_file.c_str(), std::ofstream::out | std::ofstream::binary);
     cereal::BinaryOutputArchive oarchive(rf);

     algorithms.push_back("Backbone");
     
     RankFileHeader h;
     h.set_nn( reference_genes.size() );
     h.set_ne( lt_mat.size() * ( lt_mat.size() - 1 ) / 2 );
     h.set_type( AGGR_TYPE );
     h.set_node_names( reference_genes );
     h.set_v( _XSTR(VERSION) );
     std::string cmd = argv[0];
     for(int i = 1; i < argc; i++) cmd = cmd + " " + args[i];
     h.set_cmd(cmd);
     h.set_method_names(algorithms);
     oarchive(h);

     for(uint64_t i = 1; i < lt_mat.size(); i++)
       {
         for(uint64_t j = 0; j < i; j++)
           {
             auto& v = lt_mat.at(i, j);
             if(! std::isnan(v.get_r() ) )
               {
                 auto s = v.get_r();
                 v.set_r( reverse ? max - s : s );
               }
             oarchive(v);
           }
       }
  
     rf.close();
     
    }
  catch(std::exception& except)
    {
      log(LOG_ERR) << except.what() << "\n";
      return errno;
    }
  
  return 0;

}
