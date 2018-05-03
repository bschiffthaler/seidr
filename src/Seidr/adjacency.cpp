/**
 * @file
 * @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
 * @version 0.01
 *
 * @param el_file Aggregated network in binary format."
 * @param gene_file File with gene names that the vertices can assume.
 * @return int 0 if the function succeeded, an error code otherwise.
 *
 * @section DESCRIPTION
 *
 * In order to perform some analyses with the aggregated network, it is
 * convenient to be able to transform it into a square adjacancy matrix.
 * This is the purpose of this function. The output will be useable by
 * e.g.: WGCNA.
 */

// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <cerrno>
#include <vector>
#include <string>
#include <fstream>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <stdexcept>
#include <cereal/archives/binary.hpp>

int adjacency(int argc, char * argv[]) {

  LOG_INIT_CERR();

  // Variables used by the function
  std::string el_file;
  std::string out_file;
  std::vector<std::string> genes;
  bool force = false;
  uint32_t tpos;
  uint16_t prec;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " adjacency";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Convert an aggregated network to an adjacency matrix.",
                       ' ', version);

    TCLAP::ValueArg<std::string>
    arg_outfile("o", "outfile", "Output file name", false,
                "adjacency.txt", "adjacency.txt");
    cmd.add(arg_outfile);

    TCLAP::ValueArg<uint32_t>
    arg_tindex("i", "index", "Output scores from this index", false,
               0, "last column");
    cmd.add(arg_tindex);

    TCLAP::ValueArg<uint16_t>
    arg_prec("p", "precision", "Number of significant digits to report", false,
             8, "8");
    cmd.add(arg_prec);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("infile", "Input file (aggregated gene counts)", true,
               "", "");

    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    el_file = arg_infile.getValue();
    out_file = arg_outfile.getValue();
    tpos = arg_tindex.getValue();
    prec = arg_prec.getValue();
    force = switch_force.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    std::cerr << "[Parser error]: " << except.what() << std::endl;
    return EINVAL;
  }

  try
  {
    file_can_read(el_file.c_str());

    if (file_exists(out_file) && ! force)
      throw std::runtime_error("File exists: " + out_file);

    file_can_create(out_file.c_str());
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << "[Runtime error]: " << except.what() << '\n';
    return errno;
  }

  SeidrFile inf(el_file.c_str());
  inf.open("r");

  SeidrFileHeader h;
  h.unserialize(inf);

  arma::mat res(h.attr.nodes, h.attr.nodes);
  res.fill(0);

  if (tpos == 0)
    tpos = h.attr.nalgs - 1;
  else
    tpos--;

  if (h.attr.dense)
  {
    for (uint32_t i = 1; i < h.attr.nodes; i++)
    {
      for (uint32_t j = 0; j < i; j++)
      {
        SeidrFileEdge e;
        e.unserialize(inf, h);
        if (EDGE_EXISTS(e.attr.flag))
        {
          res(i, j) = e.scores[tpos].s;
          res(j, i) = e.scores[tpos].s;
        }
      }
    }
  }
  else
  {
    for (uint64_t i = 0; i < h.attr.edges; i++)
    {
      SeidrFileEdge e;
      e.unserialize(inf, h);
      res(e.index.i, e.index.j) = e.scores[tpos].s;
      res(e.index.j, e.index.i) = e.scores[tpos].s;
    }
  }

  std::ofstream ofs(out_file.c_str(), std::ios::out);
  ofs.precision(prec);
  ofs.setf( std::ios::fixed, std:: ios::floatfield );

  for (uint32_t i = 0; i < res.n_rows; i++)
  {
    for (uint32_t j = 0; j < res.n_cols; j++)
    {
      ofs << res(i, j);
      ofs << (j == res.n_cols - 1 ? '\n' : '\t');
    }
  }

  return 0;

}
