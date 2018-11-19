// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
#include <top.h>
// External
#include <cerrno>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <queue>
#include <tclap/CmdLine.h>

typedef std::pair<SeidrFileEdge, double> expl_score_t;

// Priority queue should be sorted inverse (greater scores last) to make
// popping off the minimum element easy
class comp_greater
{
public:
  bool operator()(const expl_score_t& lhs, const expl_score_t& rhs) 
  {
    return lhs.second > rhs.second;
  }
};

// Need to derive from STL priority_queue to access the underlying vector
// for convenient iteration
class priority_queue : 
public std::priority_queue<expl_score_t, std::vector<expl_score_t>, comp_greater>
{
public:
  std::vector<expl_score_t>& vec() {return c;}
};


int top(int argc, char * argv[]) {

  logger log(std::cerr, "top");

  // Variables used by the function
  std::string el_file;
  uint32_t ntop;
  uint16_t tpos;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " top";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Extract top scoring edges", ' ', version);

    TCLAP::ValueArg<uint32_t>
    ntop_arg("n", "number", "The number of highest scoring edges to return",
             false, 1000, "1000");
    cmd.add(ntop_arg);

    TCLAP::ValueArg<uint16_t>
    arg_tindex("i", "index", "Select top edges of this index", false,
               0, "last column");
    cmd.add(arg_tindex);

    TCLAP::UnlabeledValueArg<std::string> arg_infile("in-file", "Input SeidrFile", true,
        "", "");
    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    el_file = arg_infile.getValue();
    ntop = ntop_arg.getValue();
    tpos = arg_tindex.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  try
  {
    el_file = to_canonical(el_file);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

  SeidrFile rf(el_file.c_str());
  rf.open("r");
  SeidrFileHeader h;
  h.unserialize(rf);

  rf.seek(0);

  if (tpos == 0)
    tpos = h.attr.nalgs - 1;
  else
    tpos--;

  double min = -std::numeric_limits<double>::infinity();
  double max = min;

  priority_queue pq;

  rf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    double v = e.scores[tpos].s;
    if (v > min)
    {
      pq.push(expl_score_t(e, v));
      if (v > max) max = v;
      min = pq.top().second;
    }
    if (pq.size() > ntop)
      pq.pop();
  });

  while(! pq.empty())
  {
    pq.top().first.print(h);
    pq.pop();
  }

  rf.close();

  return 0;
}
