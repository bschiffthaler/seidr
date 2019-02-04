// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
#include <describe.h>
// External
#include <cerrno>
#include <vector>
#include <string>
#include <fstream>
#include <armadillo>
#include <stdexcept>
#include <sstream>
#include <boost/program_options.hpp>
#include <iomanip>
#include <bs/describe.h>
#include <bs/histogram.h>

namespace po = boost::program_options;

void get_desc_stats(seidr_param_describe_t& params)
{
  SeidrFile sf(params.el_file.c_str());
  sf.open("r");

  SeidrFileHeader header;
  header.unserialize(sf);
  sf.seek(0);

  std::vector<BS::desc_stats> stats;
  for (uint64_t i = 0; i < header.attr.nalgs; i++)
  {
    stats.push_back(BS::desc_stats());
  }

  sf.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h)
  {
    for (uint16_t i = 0; i < h.attr.nalgs; i++)
    {
      if (std::isfinite(e.scores[i].s))
      {
        stats[i].add(e.scores[i].s);
      }
    }
  });

  sf.close();

  for (uint64_t i = 0; i < header.attr.nalgs; i++)
  {
    std::cout << i << '\t' << header.algs[i] << "\tmin\t" << stats[i].min() << '\n';
    std::cout << i << '\t' << header.algs[i] << "\tmax\t" << stats[i].max() << '\n';
    std::cout << i << '\t' << header.algs[i] << "\tmean\t" << stats[i].mean() << '\n';
    std::cout << i << '\t' << header.algs[i] << "\tn\t" << stats[i].size() << '\n';

    std::string quant;
    std::stringstream qs(params.quantiles);
    while (std::getline(qs, quant, ','))
    {
      std::cout << i << '\t' << header.algs[i] << "\tQ" << quant << '\t'
                << stats[i].quantile(std::stod(quant)) << '\n';
    }

    BS::histogram hist(stats[i], params.bins);
    std::stringstream ss;
    if (params.parseable)
    {
      hist.print_tsv(ss);
    }
    else
    {
      hist.print_vertical(ss);
    }

    std::string line;
    for (uint64_t j = 0; j < params.bins; j++)
    {
      std::getline(ss, line);
      std::cout << i << '\t' << header.algs[i] << "\tHIST\t" << line << '\n';
    }
  }

}

int describe(int argc, char * argv[]) {

  LOG_INIT_CERR();

  // Variables used by the function
  seidr_param_describe_t param;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " describe";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella;

  po::options_description opt("Options");
  opt.add_options()
  ("bins,b", po::value<uint32_t>(&param.bins)->default_value(25),
   "Number of histgram bins")
  ("parseable,p", po::bool_switch(&param.parseable)->default_value(false),
   "Print parseable histogram")
  ("quantiles,q", po::value<std::string>(&param.quantiles)->default_value("0.05,0.25,0.5,0.75,0.95"),
   "A comma seaprated string of quantiles [0.05,0.25,0.5,0.75,0.95]")
  ("help,h", "Show this help message");

  po::options_description req("Required");
  req.add_options()
  ("in-file", po::value<std::string>(&param.el_file)->required(),
   "Input SeidrFile");

  umbrella.add(opt).add(req);

  po::positional_options_description p;
  p.add("in-file", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, args).
            options(umbrella).positional(p).run(), vm);

  if (vm.count("help"))
  {
    std::cerr << umbrella << '\n';
    return EINVAL;
  }

  po::notify(vm);

  try
  {
    if (! file_can_read(param.el_file.c_str()))
      throw std::runtime_error("Cannot read input file: " + param.el_file);
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  get_desc_stats(param);

  return 0;
}