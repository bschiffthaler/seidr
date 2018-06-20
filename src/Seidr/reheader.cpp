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

  std::string tempfile = infile + ".tmp";
  SeidrFile rf(infile.c_str());
  rf.open("r");

  // Fill map with old index data
  std::map<uint32_t, std::string> old_index;
  rf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    auto hit = old_index.find(e.index.i);
    if (hit == old_index.end())
      old_index[e.index.i] = h.nodes[e.index.i];
    hit = old_index.find(e.index.j);
    if (hit == old_index.end())
      old_index[e.index.j] = h.nodes[e.index.j];
  });

  rf.close();

  // Create new vector of nodes for the header
  std::vector<std::string> new_nodes;
  for (auto& e : old_index)
    new_nodes.push_back(e.second);

  // Fill a map with new index data
  std::map<std::string, uint32_t> new_index;
  for (uint32_t i = 0; i < new_nodes.size(); i++)
    new_index[ new_nodes[i] ] = i;

  SeidrFile of(tempfile.c_str());
  of.open("w");

  SeidrFile inf(infile.c_str());
  inf.open("r");

  SeidrFileHeader h;
  h.unserialize(inf);

  auto old_nodes = h.nodes;

  // Fix nodes
  h.attr.nodes = new_nodes.size();
  h.nodes = new_nodes;
  h.serialize(of);

  inf.seek(0);
  inf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h2) {
    e.index.i = new_index[ old_index.at(e.index.i) ];
    e.index.j = new_index[ old_index.at(e.index.j) ];
    e.serialize(of, h); // Use NEW header
  });

  of.close();
  inf.close();

  rename(tempfile, infile);

  return 0;
}