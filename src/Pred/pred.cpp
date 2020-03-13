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
#include <fs.h>
#include <pred-fun.h>
#include <BSlogger.hpp>
// External
#include <string>
#include <vector>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#define MODEL_RF 0
#define MODEL_GB 1
#define MODE_TRAIN 2
#define MODE_PREDICT 3

namespace fs = boost::filesystem;

int main(int argc, char ** argv){
  
  std::string infile;
  std::string outfile;
  std::string modelfile;
  size_t classifier = MODEL_RF; //default RF
  size_t mode = MODE_PREDICT;

  LOG_INIT_CERR();
  
  // Define program options
  try 
    {
      TCLAP::CmdLine cmd("PRED C++", ' ', "0.1");

      TCLAP::ValueArg<std::string> infile_arg("i", "infile",
                                              "The data file", true,
                                              "", "string");
      cmd.add(infile_arg);

      TCLAP::ValueArg<std::string> model_arg("m", "model",
                                               "File containing ML model",
					       true, "",
                                               "string");
      cmd.add(model_arg);
      
      TCLAP::ValueArg<std::string> outfile_arg("o", "outfile",
					       "Output file path", false,
					       "", "string");
      cmd.add(outfile_arg);

      //keeping it as a string for future uses of other models
      TCLAP::ValueArg<std::string> classifier_arg("c", "classType",
					       "classifier: RF random forest, GB gradient booster"
					       , false,"RF", "string");
      cmd.add(classifier_arg);
      
      cmd.parse(argc, argv);
      
      infile = infile_arg.getValue();
      outfile = outfile_arg.getValue();
      modelfile = model_arg.getValue();

      //no outfile then infile is train and modelfile will store the model
      if(outfile == "")
	mode = MODE_TRAIN;
      
    
      if (classifier_arg.getValue()=="GB")
	classifier=MODEL_GB;
      //whatever else is RF forest for now	  
      
    } 
  catch (TCLAP::ArgException& e)  // catch any exceptions
    { 
      log(LOG_ERR) << e.error() << " for arg " << e.argId() << '\n';
      return EINVAL;
    }

  // Check all kinds of FS problems that may arise
  try
    {
      outfile = to_absolute(outfile);
      infile = to_absolute(infile);
      
      if(! file_exists(dirname(outfile)) )
        throw std::runtime_error("Directory does not exist: " +
                                 dirname(outfile));
      //TODO add modelfile things            
      if(! file_exists(infile) )
        throw std::runtime_error("File does not exist: " + infile);

      if(! regular_file(infile) )
        throw std::runtime_error("Not a regular file: " + infile);
          
      if(! file_can_read(infile) )
        throw std::runtime_error("Cannot read: " + infile);
   
    }
  catch(std::runtime_error& e)
    {
      log(LOG_ERR) << e.what() << '\n';
      return EINVAL;
    }

  try
    {
      // Get input files
      readInput(infile); //load file into global
      
      if (classifier==MODEL_RF) {
	std::cout<<"Model RF\n";
	
	if(mode == MODE_PREDICT) {
	  predictRF(modelfile, outfile);
	} else if (mode == MODE_TRAIN) {
	  trainModelRF(modelfile); //train & export model
	}
      } else if (classifier==MODEL_GB){
	std::cout<<"Model GB\n";
      
	if(mode == MODE_PREDICT) {
	  predictGB(modelfile, outfile);
	} else if (mode == MODE_TRAIN) {
	  trainModelGB(modelfile); //train & export model
	}
      }
    }
  catch(const std::runtime_error& e)
    {
      log(LOG_ERR) << e.what() << '\n';
      return 1;
    }
  catch(const std::exception& e)
    { //TODO treat other exceptions differently
      log(LOG_ERR) << e.what() << '\n';
      return 1;
    }
  return 0;

}
