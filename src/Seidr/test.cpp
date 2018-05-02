/**
 * @file
 * @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
 * @version 0.01
 *
 * @param el_file Binary representation of ranked edge list or aggregated network.
 * @param gene_file File with gene names.
 * @return int 0 if the function succeeded, an error code otherwise.
 *
 * @section DESCRIPTION
 *
 * This function simply converts the binary edge lists into human readble text.
 */
// Seidr
#include <common.h>
#include <test.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <cerrno>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <tclap/CmdLine.h>

int test(int argc, char * argv[]) {

  logger log(std::cerr, "test");

  // Variables used by the function
  std::string el_file;

  // We ignore the first argument
  const char * args[argc - 1];
  args[0] = strcat(argv[0], " test");
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Random debug functions.", ' ', version);
    TCLAP::UnlabeledValueArg<std::string> arg_infile("in-file", "Input file (aggregated gene counts)", true,
        "", "string");
    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    el_file = arg_infile.getValue();
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

  rf.each_edge([](SeidrFileEdge& e, SeidrFileHeader& h){
    e.print(h);
  });

  rf.close();

  return 0;
}
