#pragma once

#include <string>
#include <vector>
#include <set>
#include <armadillo>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <boost/archive/tmpdir.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>


void readInput(std::string filename);

void trainModelGB(std::string filename);
void predictGB(std::string modelfilename, std::string outfilename);

void trainModelRF(std::string filename);
void predictRF(std::string modelfilename, std::string outfilename);




