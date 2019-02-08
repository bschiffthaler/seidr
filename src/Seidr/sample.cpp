// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
#include <sample.h>
#include <cmath>
// External
#include <random>
#include <bs/vitter_d.h>
#include <bs/simple_sample.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

uint64_t prop_to_nedges(double prop, uint64_t n)
{
  double nd = static_cast<double>(n) * prop;
  uint64_t ret = static_cast<uint64_t>(round(nd));
  return ret;
}

void sample_wo_repl(const seidr_param_random_t& param,
                    const po::variables_map& vm)
{
  SeidrFile sf(param.el_file.c_str());
  sf.open("r");
  SeidrFileHeader header;
  header.unserialize(sf);
  sf.close();
  sf.open("r");
  uint64_t nedges;
  if (vm.count("nedges"))
  {
    if (param.nedges > header.attr.edges)
    {
      throw std::runtime_error("Can't sample " + std::to_string(param.nedges) + 
                               " edges when the network only contains " +
                               std::to_string(header.attr.edges) + " edges");
    }
    nedges = param.nedges;
  }
  else
  {
    nedges = prop_to_nedges(param.prop, header.attr.edges);
  }

  uint64_t ctr = 0;

  BS::vitter_d vd(header.attr.edges, nedges);

  uint64_t cur = vd.next();

  sf.each_edge_exit_early([&](SeidrFileEdge& e, SeidrFileHeader& h){
    if (cur == ctr)
    {
      e.print(h);
      if (vd.end())
        return true;
      else
        cur = vd.next();
    }
    ctr++;
    return false;
  });

  sf.close();
}

void sample_w_repl(const seidr_param_random_t& param,
                    const po::variables_map& vm)
{
  SeidrFile sf(param.el_file.c_str());
  sf.open("r");
  SeidrFileHeader header;
  header.unserialize(sf);
  sf.close();
  sf.open("r");
  uint64_t nedges;
  if (vm.count("nedges"))
  {
    nedges = param.nedges;
  }
  else
  {
    nedges = prop_to_nedges(param.prop, header.attr.edges);
  }

  uint64_t ctr = 0;

  BS::simple_sample<uint64_t> ss(header.attr.edges, nedges);

  uint64_t cur = ss.next();

  sf.each_edge_exit_early([&](SeidrFileEdge& e, SeidrFileHeader& h){
    while (cur == ctr)
    {
      e.print(h);
      if (ss.end())
        return true;
      else
        cur = ss.next();
    }
    ctr++;
    return false;
  });

  sf.close();
}

int random(int argc, char * argv[]) {

  logger log(std::cerr, "random");

  // Variables used by the function
  seidr_param_random_t param;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " random";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  po::options_description umbrella;

  po::options_description opt("Options");
  opt.add_options()
  ("help,h", "Show this help message")
  ("index,i", po::value<uint32_t>(&param.tpos)->default_value(0),
   "Score index to use")
  ("replacement,r", po::bool_switch(&param.replacement)->default_value(false),
   "Sample with replacement")
  ("precision,p", po::value<uint16_t>(&param.prec)->default_value(8),
   "Number of significant digits to report");

  po::options_description oneof("One of");
  oneof.add_options()
  ("nedges,n", po::value<uint64_t>(&param.nedges),
   "Number of edges to sample")
  ("fraction,f", po::value<double>(&param.prop),
   "Fraction of edges to sample");

  po::options_description req("Required");
  req.add_options()
  ("in-file", po::value<std::string>(&param.el_file)->required(),
   "Input SeidrFile");

  umbrella.add(opt).add(oneof).add(req);

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
    // Check a bunch of sanity things
    if (! file_can_read(param.el_file.c_str()))
      throw std::runtime_error("Cannot read input file: " + param.el_file);

    if (vm.count("nedges") && vm.count("fraction"))
      throw std::runtime_error("Supply either --fraction or --nedges, not both");

    if ( (! vm.count("nedges")) && (! vm.count("fraction")) )
      throw std::runtime_error("Supply either --fraction or --nedges");

    if (vm.count("fraction"))
    {
      if (!param.replacement && (param.prop < 0 || param.prop > 1))
        throw std::runtime_error("Fraction must be in [0, 1] if sampling " 
                                 "without replacement");
      if (param.replacement && param.prop < 0)
        throw std::runtime_error("Fraction must be > 0");
    }

    // Start main function
    if (! param.replacement)
    {
      sample_wo_repl(param, vm);
    }
    else
    {
      sample_w_repl(param, vm);
    }
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  return 0;
}
