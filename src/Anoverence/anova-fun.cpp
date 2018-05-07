
// Seidr
#include <anova-fun.h>
#include <BSlogger.h>
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

typedef arma::uword anovaIndex;

int geneNumber; //genes
int chipTotal; //chip;
std::vector<chip> chips; // main vector of chips rows size
std::map<std::string, gene> genes; 
std::map<int, std::vector<chip*>> experimentMap; 
std::map<int, std::vector<double>> chipValuesMap;
std::map<int, std::vector<double>> meansMap;
std::map<int, chip*> chipMap;
std::vector<std::string> geneNames; //store the gene names

/********************************************************
 FUNCTIONS
 *********************************************************/

/*
 * read expression data from file
 *
 */
void readExpressionData(std::string filename) {
  std::ifstream infile (filename);
  std::string line;
  std::getline(infile,line); 	//read first line and use it to get columns
  std::stringstream linestream(line);
  std::string data;
  while(std::getline(linestream,data,'\t')) 	// count columns
    ++geneNumber;

  //read a line and start processing, from second line
  while(std::getline(infile,line)) { 
    std::stringstream linestream(line);
    std::string data;

    while(std::getline(linestream,data,'\t')) { //each line to a vector
      chipValuesMap[chipTotal+1].push_back(std::stod(data));
    }
    chipTotal++; //line is read, add one to row count
  }
  chips.resize(chipTotal); //fixed size
}

/*
 * Read features file
 *
 */

void readFeatures(std::string filename) {
  std::ifstream infile (filename);
  std::string first_line;
  std::getline(infile, first_line); //Ignores the first line

  int chipCount=1;
  chip tempChip;
  
  for (int i=0; i<chipTotal; i++) { //for each line
    std::string currentLine;
    std::getline(infile, currentLine);
    std::vector<std::string> vLine; 

    std::stringstream ss(currentLine);
    std::string token;

    while(std::getline(ss, token, '\t')) //split string
      vLine.push_back(token);
       
    chip thisChip; //stores new chip info

    thisChip.chipNumber=chipCount;//tag the chip
    
    // 1 experiment
    thisChip.experiment= (vLine[0]=="NA") ? -1 : std::stoi(vLine[0]);

    //2&3 perturbation && perturbation levels
    if(vLine[1]!="NA") {     
      std::vector<std::string> pertTemp;
      std::vector<float> levelTemp;
      std::stringstream pertStream(vLine[1]);
      std::stringstream levelStream(vLine[2]);
      std::string piece;

      while(std::getline(pertStream,piece,',')) {
        pertTemp.push_back(piece);
        thisChip.perturbationsID.push_back(piece);
      }
      
      while(std::getline(levelStream,piece,',')) {
        if (piece=="NA") {
          for (anovaIndex pert=0; pert<pertTemp.size(); pert++)
            levelTemp.push_back(0.0); 
        } else {
          levelTemp.push_back(std::stof(piece));
        }
      }
      if(pertTemp.size()!=levelTemp.size()) {
        throw std::runtime_error("Perturbation info is missing\n");
      } else {
        for (anovaIndex i=0;i<pertTemp.size();i++) {
          sparse sTemp (pertTemp[i], levelTemp[i]);
          thisChip.perturbations.push_back(sTemp);
	      }
      }
    }
    sort(thisChip.perturbationsID.begin(), thisChip.perturbationsID.end());
    sort(thisChip.perturbations.begin(), thisChip.perturbations.end());

    // 4 treatment
    thisChip.treatment= (vLine[3]=="NA") ? -1 : std::stoi(vLine[3]);    

    //5 deletedgenes
    if(vLine[4]!="NA") { //if NA leave an empty vector
      std::stringstream delStream(vLine[4]);
      std::string delPiece;
      while(std::getline(delStream,delPiece,',')) {
        if (std::find(geneNames.begin(), geneNames.end(), delPiece)
            != geneNames.end() ){
          thisChip.deletedGenes.push_back(delPiece); 
          genes[delPiece].del++;
        }
    }
    sort(thisChip.deletedGenes.begin(), thisChip.deletedGenes.end());

    //6 overexpressed genes
    if(vLine[5]!="NA") { //if NA leave an empty vector
      std::stringstream overStream(vLine[5]);
      std::string overPiece;
      while(std::getline(overStream,overPiece,',')) {
        if (std::find(geneNames.begin(), geneNames.end(), delPiece)
            != geneNames.end() ) {
          thisChip.overexpressedGenes.push_back(overPiece);
          genes[overPiece].over++;
        }
      }
    }
    sort(thisChip.overexpressedGenes.begin(),
         thisChip.overexpressedGenes.end());

    // 7 time
    thisChip.time= (vLine[6]=="NA") ? 0 : std::stof(vLine[6]);

    // 8 time
    thisChip.replicaNumber= (vLine[7]=="NA") ? 0 : std::stoi(vLine[7]); //NEW

    //store this chip info in temp waiting for the data from the replicates
    chips[i]=thisChip; 
    chipCount++;
  } //end for line
  std::sort(chips.begin(), chips.end(), compareChips); //sort the chips vector
} //end function


/*
 * read gene files for gene name
 *
 */
void readGeneFile(std::string filename) {
  std::ifstream infile (filename);
  std::string line;

  while(std::getline(infile, line)) 
    geneNames.push_back(line);

}


/*
 * Creates map of pointers for better performance
 *
 */
void mapChipVector() {

  chip* tempChip{};
  
  for (auto & c:chips){ //they are ordered
    chipMap[c.chipNumber]=&c;
    if(c.replicaNumber==1) {
      tempChip=&c;
      experimentMap[c.experiment].push_back(&c);
    } else {
      tempChip->totalRepeats++;
    }
  }
}

/*
 * Calculate the mean
 */

void calculateMean(){

  int repeats;
  int kIndex;	
  double meanTemp;

  for (anovaIndex k=0; k<chips.size(); k++) {
    if (chips[k].replicaNumber!=1) //don't look at the replicas
      continue;

    kIndex=chips[k].chipNumber;	
    repeats=chips[k].totalRepeats;

    for (int g=0; g<geneNumber; g++){
      meanTemp=0.0;
      for (int i=0; i<repeats; i++) {
        meanTemp += chipValuesMap[kIndex+i][g];
      }
      meansMap[kIndex].push_back(meanTemp/repeats);
    } //end for each gene
  } //end for each chip
}

/*
 * Selects a control condition c for a measurement k 
 * following requirements 1a, 1b or 2
 * PRECONDITION - chips vector MUST be ordered
 */
void selectPairs() {

  std::vector<chip*>controlsTemp; //keep for the requirement 2

  for (auto &itr : experimentMap ) { //for each pair
    for(auto &k :itr.second) {
      controlsTemp.clear();
      bool equalFound=false;
      
      for(auto &c :itr.second) { //TODO maybe we can top the iteration until k index
        if (k == c) { //don't compare itself
          continue;
        }
	
        if (totalPerturbed(c,k )>0)  //k must not have a subset of c
          continue;
        if(isPerturbationEqualSet(c,k) && c->time < k->time) {
          if (!equalFound) {
            k->controlsPtr.push_back(c);
            equalFound=true;
            continue;
          }
        }
        if (totalPerturbed(k,c)>0) { //cis a subset of k
          controlsTemp.push_back(c);
        }
      }//end for each control
      
      chip* bestChip;

      for (anovaIndex can =0; can<controlsTemp.size(); can++) {
        if (controlsTemp.size()==1) { //only one candidate? then is the best
          k->controlsPtr.push_back(controlsTemp[can]);
          break;
        }

        if (can==0){ //skip the 1st one (can't be compared)
          bestChip=controlsTemp[can];
          continue;
        }
        //this an previous are from the same group
        if(isPerturbationEqualSet(controlsTemp[can-1], controlsTemp[can])) {
          if ( std::fabs(k->time - controlsTemp[can]->time)
               < std::fabs(k->time - bestChip->time) ) {
	          bestChip=controlsTemp[can];	    
	        }
	      } else {
	        k->controlsPtr.push_back(bestChip);
	        bestChip = controlsTemp[can];
	      } 

        //whatever happens the last must be pushed
        if (can==controlsTemp.size()-1) {
          k->controlsPtr.push_back(bestChip);
        }
	
	
      } //end for
    } //end k for
  } //end experiment group for
} //end function



double twoWayAnova (int genePair[2], float wko) {

  double n=0; //N
  double wsum=0; //wsum
  double Xj=0; //Xj
  double Xijk=0; //Xijk
  double wgt, sum, Fijk, Pijk;
  double x, SSa,SSt;


  for (anovaIndex k=0; k<chips.size(); ++k) {
    if (chips[k].replicaNumber!=1) {//don't look at the replicas
      continue;
    }
    int repeats=chips[k].totalRepeats; //for performance
    for (auto &control : chips[k].controlsPtr ) { //for each control of that k 
      wgt=0;
      sum=0;
      

      for(int j=0; j<2; j++) { //for each gene
        double mean = meansMap[control->chipNumber][genePair[j]];
        for (int i=0; i<repeats;++i) {
          int kIndex=chips[k+i].chipNumber;
          Fijk = chipValuesMap[kIndex][genePair[j]] - mean;
          Pijk = (j==1 && isPerturbed(genePair[j], &chips[k]) 
                 && !isPerturbed(genePair[j], control)) ? wko: 1;

          wgt += Pijk;
          sum += Pijk*Fijk;
          Xijk += Pijk*Fijk*Fijk;
        }// end for each replica
      }//end for each gene

      n+=wgt;
      wsum+=sum;
      Xj+=sum*sum/wgt;
    } //end for each control
  } //end for each condition
   
  x = wsum*wsum / n;
  SSa = Xj - x;
  SSt = Xijk - x;
  return SSa / SSt;
}

void analizeGenePairs(std::string filename, float wko){
  LOG_INIT_CERR();
  log(LOG_INFO) << "Analyzing gene pairs...\n";
  
  int endColumn=1; //keeps track for the diagonal matrix
  int genePair[2]={0,0};
  std::ofstream outfile;
  outfile.open(filename);
  
  int elements=0;

  for (int i=1; i<geneNumber; ++i){ //skip the first line
    genePair[0]=i; //rows starging at gene 2
    for (int j=0; j<endColumn; ++j){
      genePair[1]=j; //columns starting at gene 1
      if (j == geneNumber)
        break;
      outfile<<twoWayAnova(genePair,wko);
      if (j < endColumn-1)
        outfile<<"\t";
      elements++;
    }

    endColumn++;
    outfile<<"\n";
  }
  outfile.close();
  
  log(LOG_INFO) << "Completed\n";
}

void analizePartialPairs(std::vector<std::string> targets,
       std::string filename, float wko) {
  LOG_INIT_CERR();
  log(LOG_INFO) << "Analyzing gene pairs...\n";

  int genePair[2]={0,0};
  std::ofstream outfile;
  outfile.open(filename);

  
  for (int j=0; j<geneNumber; ++j) {
    genePair[0]=j;
    for (auto t:targets) {
      auto pos = std::distance(geneNames.begin(), find(geneNames.begin(),
                 geneNames.end(), t));
      genePair[1]=pos;
      if(genePair[0]==genePair[1])
         continue;
      outfile<<geneNames[genePair[1]]<<"\t"<<geneNames[genePair[0]]
	     <<"\t"<<twoWayAnova(genePair, wko)<<"\n";

    }   
  }

  outfile.close();
  log(LOG_INFO) << "Completed\n";
}


/****************************************************
 * AUXILIAR FUNCTIONS
 ***************************************************/

/*
 * countPerturbed
 * return the number of exclusive conditions in first set
 */


int totalPerturbed(chip* c1, chip* c2){

  return
    countPerturbed(c1->perturbationsID, c2->perturbationsID) +
    countPerturbed(c1->deletedGenes, c2->deletedGenes) +
    countPerturbed(c1->overexpressedGenes, c2->overexpressedGenes) +
    (c1->treatment>=0 && c2->treatment<0);
}

/*
 * isPerturbationEqual
 * returns true if chip1 has additional conditions
 * perturbations,deletions, overexpressions or different treatments
 * if conditions are the same returns false
 *
 */
bool isPerturbationEqualSet (chip* c1, chip* c2) {

  return
    (c2->perturbationsID == c1->perturbationsID) &&
    (c2->deletedGenes == c1->deletedGenes) &&
    (c2->overexpressedGenes == c1->overexpressedGenes) &&
    (c1->treatment == c2->treatment);
}

/*
 * isPerturbed
 * returns true if a gene is deleted or overexpressed in a chip
 */
bool isPerturbed(int gene, chip* c){


  for (std::string dGene: c->deletedGenes){
    //std::cout<<"deleted "<<dGene<<"\n";
    if (dGene==geneNames[gene]) {
      return true;
    }
  }
  for (std::string oGene: c->overexpressedGenes){
    //std::cout<<"deleted "<<oGene<<"\n";
    if (oGene==geneNames[gene]) {
      return true;
    }
  }
  return false;
}

/*
 * compareChips
 * returns true chip1 should be before chip2
 * Used with sort function
 */
bool compareChips(chip c1, chip c2) {
  LOG_INIT_CERR();
  // Crit 1
  if (c1.experiment < c2.experiment) //less experiment first
    return true;
  else if (c1.experiment > c2.experiment) //less experiment first
    return false;

  // Crit 2 - perturbations
  if(c1.perturbations.size() < c2.perturbations.size())
    return true;
  else if (c1.perturbations.size() > c2.perturbations.size())
    return false;
  else if (c1.perturbations.size() == c2.perturbations.size()) {
    if(c1.perturbations.size() != 0) {
      for (anovaIndex i=0; i<c1.perturbations.size(); i++) {
        //compare each pair
        if (c1.perturbations[i].index == c2.perturbations[i].index) {
          // the lowest level wins
          if (c1.perturbations[i].value < c2.perturbations[i].value) { 
            return true;
          } else if (c1.perturbations[i].value == c2.perturbations[i].value) {
            continue;
          } else {
            return false;
          }
        } else if (c1.perturbations[i].index < c2.perturbations[i].index) {
          return true;
        } else {
          return false;
        }
      }//end for
    } 
  } //end of criteria 2

  //Crit 3
  if(c1.treatment < c2.treatment)
    return true;
  else if (c1.treatment > c2.treatment)
    return false;

  // Crit 4
  if(c1.deletedGenes.size() < c2.deletedGenes.size())
    return true;
  else if (c1.deletedGenes.size() > c2.deletedGenes.size())
    return false;
  else if (c1.deletedGenes.size() == c2.deletedGenes.size()) {
    if(c1.deletedGenes.size() != 0) {
      for (anovaIndex i=0; i<c1.deletedGenes.size(); i++) {	//compare each pair
        if (c1.deletedGenes[i] == c2.deletedGenes[i])
          continue; //keep searching
        else if (c1.deletedGenes[i] < c2.deletedGenes[i])
          return true;
        else 
          return false;
      }
    }
  }

  // Crit 5
  if(c1.overexpressedGenes.size() < c2.overexpressedGenes.size())
    return true;
  else if (c1.overexpressedGenes.size() > c2.overexpressedGenes.size())
    return false;
  else if (c1.overexpressedGenes.size() == c2.overexpressedGenes.size()) {
    if(c1.overexpressedGenes.size() != 0) {
      for (anovaIndex i=0; i<c1.overexpressedGenes.size(); i++) {
	if (c1.overexpressedGenes[i] == c2.overexpressedGenes[i])
	  continue; //keep searching
	else if (c1.overexpressedGenes[i] < c2.overexpressedGenes[i])
	  return true;
	else return false;
      }
    }
  }

  //Crit 6
  if(c1.time < c2.time)
    return true;
  else if (c1.time > c2.time)
    return false;

  //Crit 7
  if(c1.replicaNumber < c2.replicaNumber)
    return true;
  else if (c1.replicaNumber > c2.replicaNumber)
    return false;

  //if it gets here both chips are the same, this should not happen...
  log(LOG_WARN)<< "Chip "
           << c1.chipNumber
           << " and chip "
           << c2.chipNumber
           << " are repeated.\n";
  return false;
}
