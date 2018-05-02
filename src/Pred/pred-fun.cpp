
// Seidr
#include <pred-fun.h>
#include <BSlogger.h>
#include <common.h>

// External
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <typeinfo>
#include <map>
#include <set>
#include <armadillo>
#include <Classifier.h>
#include "ForestProbability.h"
#include "Data.h"
#include "DataDouble.h"

// Gradient Booster parameters
seidr_uword_t gbEstimators = 2000;
seidr_uword_t gbMaxDepth = 5;

// Random Forest parameters
seidr_uword_t rfEstimators = 91; //def 20
seidr_uword_t min_node_size = 1;
double alpha = 1.0;
double minprop = 1.0;

arma::mat inputMatrix; //matrix to hold data

typedef std::vector<float> stdFloatVec; //needed for arma::conv_to

/*
 * Class PredForest
 *
 */

class PredForest : public ForestProbability {
public:
  PredForest(std::string dependent_variable_name, Data* input_data,
	     size_t ntree, size_t mtry, size_t min_node_size, double alpha,
	     double minprop, std::ostream& logstream, std::string forestfile,
	     bool test);
};

//constructor to load from file or save model
PredForest::PredForest(std::string dependent_variable_name, Data* input_data,
		       size_t ntree, size_t mtry, size_t min_node_size,
		       double alpha, double minprop, std::ostream& logstream,
		       std::string forestfile, bool test){

  std::vector<std::string> catvars;
  MemoryMode M = MEM_DOUBLE;
  ImportanceMode I = IMP_GINI;
  SplitRule S = LOGRANK;
  PredictionType P = RESPONSE;
  
  init(dependent_variable_name, M, input_data, mtry,
       //output_prefix, num_tree, seed, threads, importance mode
       forestfile, ntree, 314159265, 1, I,
       //min node size (0=auto), status variable, prediction mode,
       // sample with replacement
       min_node_size, "", test, true,
       //categorical vars, save memory, split rule (1 = default/variance)
       catvars, false, S,
       //predict all, sample fraction, alpha, min prop, holdout
       false, 1, alpha, minprop, false,
       //prediciton mode (1=response), number of random splits
       P, 1);
  
  verbose_out = &logstream; 
  
  if (test) {
    loadFromFile(forestfile);
  }
 
}

/*
 *********************************************************************
 * FUNCTIONS
 *********************************************************************
 */

/**
 * Reads data from file and put it in armadillo matrix
 * @param filename The name of the file to be processed
 */
void readInput(std::string filename) {
  //read file as arma matrix, file must be just numbers, no header.
  inputMatrix.load(filename);
}


/**
 * Train the model using Gradient Booster classifier
 * @param filename The name of the file to use as model
 */
void trainModelGB(std::string filename) {
  std::cout<<"Train model: "<<filename<<"\n";
  FastBDT::Classifier classifier;
  seidr_uword_t lastColumn = inputMatrix.n_cols-1; //for eficiency 
  std::vector<bool> classType(inputMatrix.n_rows); //store the class
  std::vector<float> classWeight(inputMatrix.n_rows,1.0); //store the weight
  std::vector<std::vector<float>> dataVector(lastColumn); //transpose matrix

  //TODO - calculate weigths and populate classWeight when needed
  
  // read the class column and turn it into bool vector
  for (seidr_uword_t i = 0; i<(inputMatrix.n_rows);++i) {
    classType[i]=(inputMatrix.at(i,lastColumn)==0)?false: true;
  }
  
  //convert matrix to vector of vector (each vector is a column!!)
  for (seidr_uword_t i=0; i<lastColumn;++i){
    dataVector[i]= arma::conv_to<stdFloatVec>::from(inputMatrix.col(i));
  }

  //no other parameters need tunning at the moment
  classifier.SetBinning({5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5}); // TODO, get more info
  classifier.SetNTrees(gbEstimators); // default is 100
  classifier.SetDepth(gbMaxDepth); // default is 3 
  classifier.SetTransform2Probability(true); //gives back probability, not class

  //run the classifier
  classifier.fit(dataVector, classType, classWeight);

  std::fstream out_stream(filename, std::ios_base::out | std::ios_base::trunc);
  out_stream << classifier << std::endl;
  out_stream.close();
}

/**
 * Train the model using Gradient Booster classifier
 * @param modelfilename The name of the GB model
 * @param outfilename The name of the file that will store output
 */
void predictGB(std::string modelfilename, std::string outfilename) {
  std::cout<<"Predict to: "<<outfilename<<"\n";
  std::fstream inStream(modelfilename, std::ios_base::in);
  FastBDT::Classifier classifier(inStream);  
  std::ofstream outfile;

  outfile.open(outfilename);
 
  std::vector<float> vectorTemp;
  for (seidr_uword_t i=0; i<inputMatrix.n_rows-1; ++i) {
    vectorTemp = arma::conv_to<stdFloatVec>::from(inputMatrix.row(i));
    outfile<<classifier.predict(vectorTemp)<<"\n"; //predict each sample
  }
  outfile.close();
}

/**
 * Train the a model using a Random Forest classifier
 * @param filename The name of the file to be processed
 */
void trainModelRF(std::string filename) {
  std::cout<<"Train model RF: "<<filename<<"\n";
  seidr_uword_t mtry = sqrt(inputMatrix.n_cols); //sqrt of variables
  std::ostream nullstream(0); //for ranger compatibility
  // a string vector is needed for the DataDouble object, it would be filled
  // with as many varX string as the number of variables -1, the last column
  // is the class type not a variable.
  std::vector<std::string> variableNames; 
  for (seidr_uword_t i = 0; i<(inputMatrix.n_cols-1);++i) {
    std::string str("var" + std::to_string(i));
    variableNames.push_back(str);
  }
  variableNames.push_back("ClassType"); //always the last column

  Data *data = new DataDouble(inputMatrix.memptr(), variableNames, inputMatrix.n_rows, inputMatrix.n_cols);
    
  PredForest forest("ClassType", data, rfEstimators, mtry, min_node_size,
		    alpha, minprop, nullstream, filename,false); //create PredForest object
    
  forest.run(false); //generate the forest

  forest.saveToFile(); //export the model to a file
  //the default name has .forest at the end, can't be changed.
    
}


/*
 * Train the model using Random Forest classifier
 * @param modelfilename The name of the GB model
 * @param outfilename The name of the file that will store output
 */
void predictRF(std::string modelfilename, std::string outfilename) {
  std::cout<<"Predicting to RF: "<<outfilename<<"\n";

  //TODO join this chunk of code with the trainRF() until forest.run(false)
  // keep it like this at the moment
  std::ostream nullstream(0); //for ranger compatibility
  seidr_uword_t mtry = sqrt(inputMatrix.n_cols); //sqrt of features
  modelfilename.append(".forest"); //append the extension
  std::ofstream outfile;
  
  outfile.open(outfilename);

  std::vector<std::string> variableNames; 
  for (seidr_uword_t i = 0; i<(inputMatrix.n_cols-1);++i) {
    std::string str("var" + std::to_string(i));
    variableNames.push_back(str);
  }
  variableNames.push_back("ClassType"); //always the last column
  
  Data* data = new DataDouble(inputMatrix.memptr(), variableNames,
			inputMatrix.n_rows, inputMatrix.n_cols);
    
  PredForest forest("ClassType", data, rfEstimators, mtry, min_node_size,
		      alpha, minprop, nullstream,modelfilename, true);

  forest.run(false); //generate the forest

  //predictions are stored in a 3D estructure
  std::vector<std::vector<std::vector<double>>> predictions = forest.getPredictions();
  for (seidr_uword_t i = 0; i < predictions[0].size(); ++i) {
    outfile << predictions[0][i][1] << "\n";
  }
    
  outfile.close();
}

