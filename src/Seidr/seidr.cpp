// Seidr
#include <common.h>
#include <import.h>
#include <viewRanks.h>
#include <aggregate.h>
#include <adjacency.h>
#include <backbone.h>
#include <index.h>
#include <threshold.h>
#include <roc.h>
#include <convert.h>
#include <backbone.h>
#include <compare_clusters.h>
#include <compare.h>
#include <resolve.h>
#include <reheader.h>
#include <stats.h>
#include <BSlogger.h>
#include <neighbours.h>
#include <test.h>
#include <top.h>
// External
#include <string>
#include <vector>
#include <stdexcept>

std::string usage_msg =
  "\n"
  "Gene network utility for format conversion, ranking and aggregation.\n"
  "--------------------------------------------------------------------\n"
  ""
  "seidr adjacency           \t Convert aggregated network to adjacency\n"
  "                          \t matrix.\n"
  "seidr aggregate           \t Aggregate a set of SeidrFiles into a crowd\n"
  "                          \t network.\n"
  "seidr cluster-enrichment  \t Test wether members of clusters in two\n"
  "                          \t networks overlap significantly or extract\n"
  "                          \t clusters.\n"
  "seidr compare             \t Compare two networks for shared/unique\n"
  "                          \t edges.\n"
  "seidr convert             \t Interconvert various text based formats.\n"
  "seidr import              \t Convert network text files to SeidrFiles.\n"
  "seidr index               \t Create index for SeidrFiles.\n"
  "seidr neighbours          \t Extract N neighbours of all nodes or a list\n"
  "                          \t of nodes in a SeidrFile.\n"
  "seidr reheader            \t Modify SeidrFile headers.\n"
  "seidr resolve             \t Resolve node indices in text file to node\n"
  "                          \t names.\n"
  "seidr roc                 \t Compute ROC curves of predictions in \n"
  "                          \t SeidrFiles given true edges.\n"
  "seidr stats               \t Compute a variety of node and edge centrality\n"
  "                          \t measure of a SeidrFile.\n"
  "seidr threshold           \t Calculate network threshold based on scale\n"
  "                          \t free fit and transitivity.\n"
  "seidr view                \t View or filter SeidrFiles.\n\n"
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
    else if (task == "stats")
    {
      ret = stats(argc, argv);
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
    return errno;
  }
  return ret;
}
