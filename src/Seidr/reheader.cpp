#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <reheader.h>
#include <BSlogger.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <tclap/CmdLine.h>

int reheader(int argc, char ** argv)
{
  logger log(std::cerr, "reheader");

  std::string infile;
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " stats";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Modify SeidrFile headers.\n"
                       "Currently only drops disconnected nodes and resets\n"
                       "stats.", ' ', version);

    TCLAP::UnlabeledValueArg<std::string> 
    arg_infile("in-file", "Input file", true,
        "", "string");

    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  std::map<uint32_t, uint32_t> reindex;
  std::set<std::string> subnodes;
  std::vector<std::string> newnodes;

  std::string tempfile = infile + ".tmp";

  SeidrFile rf(infile.c_str());
  rf.open("r");

  rf.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h){
    reindex[e.index.i] = 0;
    reindex[e.index.j] = 0;
    subnodes.insert(h.nodes[e.index.i]);
    subnodes.insert(h.nodes[e.index.j]);
  });

  uint32_t ctr = 0;
  for(auto& i : reindex)
    i.second = ctr++;

  for(auto& i : subnodes)
    newnodes.push_back(i);

  rf.close();

  SeidrFile of(tempfile.c_str());
  of.open("w");

  SeidrFile inf(infile.c_str());
  inf.open("r");

  SeidrFileHeader h;
  h.unserialize(inf);

  // Fix nodes
  h.attr.nodes = newnodes.size();
  h.nodes = newnodes;
  h.serialize(of);

  inf.seek(0);
  inf.each_edge([&](SeidrFileEdge& e, SeidrFileHeader& h2){
    e.index.i = reindex[e.index.i];
    e.index.j = reindex[e.index.j];
    e.serialize(of, h); // Use NEW header
  });

  of.close();
  inf.close();

  rename(tempfile, infile);

  return 0;
}