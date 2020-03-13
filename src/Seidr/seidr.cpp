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
#include <import.h>
#include <viewRanks.h>
#include <aggregate.h>
#include <adjacency.h>
#include <asp.h>
#include <backbone.h>
#include <describe.h>
#include <index.h>
#include <graphstats.h>
#include <threshold.h>
#include <roc.h>
#include <convert.h>
#include <compare_clusters.h>
#include <compare.h>
#include <resolve.h>
#include <reheader.h>
#include <stats.h>
#include <tau.h>
#include <BSlogger.hpp>
#include <neighbours.h>
#include <sample.h>
#include <test.h>
#include <top.h>
// External
#include <string>
#include <vector>
#include <stdexcept>

std::string usage_msg =
  "Gene network utility for format conversion, ranking and aggregation.\n\n"
  "[Create crowd networks]\n"
  "  import                    \t Convert network text files to SeidrFiles.\n"
  "  aggregate                 \t Aggregate a set of SeidrFiles into a crowd\n"
  "                            \t network.\n"
  "\n"
  "[Filter, threshold or search SeidrFiles]\n"
  "  backbone                  \t Calculate network backbone and filter edges\n"
  "                            \t based on noise corrected backbone measure.\n"
  "  index                     \t Create index for SeidrFiles.\n"
  "  neighbours                \t Extract N first degree neighbours of all nodes\n"
  "                            \t or a list of nodes in a SeidrFile.\n"
  "  sample                    \t Sample random edges from a SeidrFile.\n"
  "  threshold                 \t Calculate network threshold based on scale\n"
  "                            \t free fit and transitivity.\n"
  "  view                      \t View, filter or search SeidrFiles.\n\n"
  "\n"
  "[Calculate network statistics]\n"
  "  stats                     \t Compute node and edge centrality\n"
  "  graphstats                \t Calculate summary statistics of the network\n"
  "\n"
  "[Format conversion]\n"
  "  adjacency                 \t Convert a SeidrFile to an adjacency\n"
  "                            \t matrix.\n"
  "  convert                   \t Interconvert various text based formats.\n"
  "  resolve                   \t Convert node indices in text file to node\n"
  "                            \t names.\n"
  "\n"
  "[Compare networks]\n"
  "  cluster-enrichment        \t Test wether members of clusters in two\n"
  "                            \t networks overlap significantly or extract\n"
  "                            \t clusters.\n"
  "  compare                   \t Compare two networks for shared/unique\n"
  "                            \t edges.\n"
  "\n"
  "[Evaluate networks]\n"
  "  roc                       \t Compute ROC curves of predictions in \n"
  "                            \t SeidrFiles given true edges.\n"
  "\n"
  "[Other utility]\n"
  "  reheader                  \t Modify SeidrFile headers.\n"
  "\n"
  "Version " + version + "\n";

int main(int argc, char* argv[])
{
  logger log(std::cerr, "seidr");

  try
  {
    if (argc < 2)
    {
      std::cerr << usage_msg << std::endl;
      throw std::invalid_argument("Too few arguments.");
    }
  }
  catch (std::invalid_argument& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return EPERM;
  }

  std::string task = argv[1];

  int ret;

  try
  {
    if (task == "adjacency")
    {
      ret = adjacency(argc, argv);
    }
    else if (task == "aggregate")
    {
      ret = aggregate(argc, argv);
    }
    else if (task == "asp")
    {
      ret = asp(argc, argv);
    }
    else if (task == "backbone")
    {
      ret = backbone(argc, argv);
    }
    else if (task == "cluster-enrichment")
    {
      ret = cluster_enrichment(argc, argv);
    }
    else if (task == "compare")
    {
      ret = compare(argc, argv);
    }
    else if (task == "convert")
    {
      ret = convert(argc, argv);
    }
    else if (task == "describe")
    {
      ret = describe(argc, argv);
    }
    else if (task == "graphstats")
    {
      ret = graphstats(argc, argv);
    }
    else if (task == "import")
    {
      ret = import(argc, argv);
    }
    else if (task == "index")
    {
      ret = index(argc, argv);
    }
    else if (task == "neighbours")
    {
      ret = neighbours(argc, argv);
    }
    else if (task == "reheader")
    {
      ret = reheader(argc, argv);
    }
    else if (task == "resolve")
    {
      ret = resolve(argc, argv);
    }   
    else if (task == "roc")
    {
      ret = roc(argc, argv);
    } 
    else if (task == "sample")
    {
      ret = sample(argc, argv);
    }
    else if (task == "stats")
    {
      ret = stats(argc, argv);
    }
    else if (task == "tau")
    {
      ret = tau(argc, argv);
    }
    else if (task == "threshold")
    {
      ret = threshold(argc, argv);
    }
    else if (task == "top")
    {
      ret = top(argc, argv);
    }
    else if (task == "view")
    {
      ret = view(argc, argv);
    }    
    else if (task == "test")
    {
      ret = test(argc, argv);
    }
    
    else
    {
      log(LOG_ERR) << usage_msg << '\n';
      throw std::invalid_argument("Unrecognized task.");
    }
  }
  catch (std::exception& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return 1;
  }
  return ret;
}
