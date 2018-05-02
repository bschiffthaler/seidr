// Seidr
#include <common.h>
#include <resolve.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <cerrno>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <tclap/CmdLine.h>

void resolve_im(SeidrFileHeader& h, std::string& fname)
{

  std::ifstream ifs(fname.c_str(), std::ios::in);
  std::string line;
  uint32_t sc_max = 0;
  while (std::getline(ifs, line))
  {
    if (line.at(0) == '#') continue;
    std::string field;
    std::stringstream ss(line);
    ss >> field;
    uint32_t sc = 0;
    for (char& c : field)
      if (c == ':')
        sc++;
    sc_max = sc > sc_max ? sc : sc_max;
  }

  ifs.clear();
  ifs.seekg(0, std::ios_base::beg);

  while (std::getline(ifs, line))
  {
    if (line.at(0) == '#') continue;
    std::string field;
    std::stringstream ss(line);
    uint32_t ctr = 0;
    while (std::getline(ss, field, ' '))
    {
      if (ctr == 0)
      {
        std::cout << field << '\t';
        uint32_t sc = 0;
        for (char c : field)
          if (c == ':')
          {
            std::cout << '\t';
            sc++;
          }
          else
            std::cout << c;
        while (sc < sc_max)
        {
          std::cout << '\t' << "NA";
          sc++;
        }
        std::cout << '\t';
        ctr++;
      }
      else if (ctr == 3)
      {
        auto index = std::stoul(field);
        std::cout << h.nodes[index] << '\n';
        ctr = 0;
      }
      else
      {
        std::cout << field << '\t';
        ctr++;
      }
    }
  }
}

int resolve(int argc, char * argv[])
{

  logger log(std::cerr, "resolve");

  // Variables used by the function
  std::string file_in;
  std::string format;
  std::string sf;

  // We ignore the first argument, the function name
  const char * args[argc - 1];
  std::string a0(argv[0]);
  a0 += " aggregate";
  const char * a1 = a0.c_str();
  args[0] = a1;
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  std::vector<std::string> format_constraints{"infomap"};

  try
  {
    TCLAP::ValuesConstraint<std::string> constraints_m(format_constraints);
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Resolve node indices in text file to node names.",
                       ' ', version);

    TCLAP::ValueArg<std::string>
    arg_format("f", "format", "File format to resolve", false,
               "infomap", &constraints_m);
    cmd.add(arg_format);

    TCLAP::ValueArg<std::string>
    arg_seidrfile("s", "seidr-file",
                  "Seidr file which should be used to resolve input", false,
                  "", "<file>");
    cmd.add(arg_seidrfile);

    TCLAP::UnlabeledValueArg<std::string>
    arg_input_file("infiles", "Input file",
                   true, "", "string");
    cmd.add(arg_input_file);

    // Parse arguments
    cmd.parse(argc, args);
    format = arg_format.getValue();
    file_in = arg_input_file.getValue();
    sf = arg_seidrfile.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return EINVAL;
  }

  if (format == "infomap")
  {
    SeidrFile fin(sf.c_str());
    fin.open("r");
    SeidrFileHeader h;
    h.unserialize(fin);
    resolve_im(h, file_in);
  }

  return 0;
}
