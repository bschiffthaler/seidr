#include <common.h>
#include <index.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include <cmath>

int index(int argc, char ** argv)
{
  logger log(std::cerr, "index");
  std::string infile;
  bool force = false;

  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " index";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Create index for SeidrFiles", ' ', version);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("in-file", "Input file", true, "", "string");

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    infile = arg_infile.getValue();
    force = switch_force.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  std::string outfile = infile + ".sfi";

  try
  {
    file_can_read(infile.c_str());
    infile = to_canonical(infile);
    if (file_exists(outfile) && ! force)
      throw std::runtime_error("File exists: " + outfile);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

  SeidrFile inf(infile.c_str());

  log(LOG_INFO) << "Building index\n";

  SeidrFileIndex sfi;
  sfi.build(infile);

  log(LOG_INFO) << "Writing index to " << outfile << '\n';

  SeidrFile otf(outfile.c_str());
  otf.open("w");
  sfi.serialize(otf);
  otf.close();
  return 0;
}