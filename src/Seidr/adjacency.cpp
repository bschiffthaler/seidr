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
 * convenient to be able to transform it into a square adjacency matrix.
 * This is the purpose of this function. The output will be useable by
 * e.g.: WGCNA.
 */

// Seidr
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
#include <adjacency.h>
// External
#include <cerrno>
#include <vector>
#include <string>
#include <fstream>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <stdexcept>

void make_full(const seidr_param_adjacency_t& param)
{
  SeidrFile sf(param.el_file.c_str());
  sf.open("r");
  SeidrFileHeader h;
  h.unserialize(sf);
  uint32_t si = param.tpos;
  if (si == 0)
    si = h.attr.nalgs - 1;
  else
    si -= 1;
  std::shared_ptr<std::ostream> out;
  if (param.out_file == "-")
    out = std::shared_ptr < std::ostream > (&std::cout, [](void*) {});
  else
    out = std::shared_ptr<std::ostream>
          (new std::ofstream(param.out_file.c_str()));

  (*out) << std::setprecision(param.prec);

  arma::mat res(h.attr.nodes, h.attr.nodes);
  res.fill(std::numeric_limits<double>::quiet_NaN());
  if (h.attr.dense)
  {
    for (uint32_t i = 1; i < h.attr.nodes; i++)
    {
      for (uint32_t j = 0; j < i; j++)
      {
        SeidrFileEdge e;
        e.unserialize(sf, h);
        if (EDGE_EXISTS(e.attr.flag))
        {
          res(i, j) = e.scores[si].s;
          res(j, i) = e.scores[si].s;
        }
      }
    }
  }
  else
  {
    for (uint64_t i = 0; i < h.attr.edges; i++)
    {
      SeidrFileEdge e;
      e.unserialize(sf, h);
      res(e.index.i, e.index.j) = e.scores[si].s;
      res(e.index.j, e.index.i) = e.scores[si].s;
    }
  }

  for (uint32_t i = 0; i < res.n_rows; i++)
  {
    for (uint32_t j = 0; j < res.n_cols; j++)
    {
      if (std::isfinite(res(i, j)))
        *out << res(i, j);
      else
        *out << param.fill;
      *out << (j == res.n_cols - 1 ? '\n' : '\t');
    }
  }
  sf.close();
}

void make_lower(const seidr_param_adjacency_t& param)
{
  SeidrFile sf(param.el_file.c_str());
  sf.open("r");
  SeidrFileHeader h;
  h.unserialize(sf);
  uint32_t si = param.tpos;
  if (si == 0)
    si = h.attr.nalgs - 1;
  else
    si -= 1;
  std::shared_ptr<std::ostream> out;
  if (param.out_file == "-")
    out = std::shared_ptr < std::ostream > (&std::cout, [](void*) {});
  else
    out = std::shared_ptr<std::ostream>
          (new std::ofstream(param.out_file.c_str()));

  (*out) << std::setprecision(param.prec);

  if (h.attr.dense)
  {
    if (param.diag)
    {
      for (uint32_t i = 0; i < h.attr.nodes; i++)
      {
        for (uint32_t j = 0; j <= i; j++)
        {
          if (j == i)
          {
            (*out) << param.fill << '\n';
          }
          else
          {
            SeidrFileEdge e;
            e.unserialize(sf, h);
            if (EDGE_EXISTS(e.attr.flag))
            {
              (*out) << e.scores[si].s << '\t';
            }
            else
            {
              (*out) << param.fill << '\t';
            }
          }
        }
      }
    }
    else // no diagonal
    {
      for (uint32_t i = 1; i < h.attr.nodes; i++)
      {
        for (uint32_t j = 0; j < i; j++)
        {
          SeidrFileEdge e;
          e.unserialize(sf, h);
          if (EDGE_EXISTS(e.attr.flag))
          {
            (*out) << e.scores[si].s;
          }
          else
          {
            (*out) << param.fill;
          }
          *out << (j == i - 1 ? '\n' : '\t');
        }
      }
    }
  }
  else //sparse
  {
    if (param.diag)
    {
      SeidrFileEdge e;
      e.unserialize(sf, h);
      uint64_t ctr = 1;
      for (uint32_t i = 0; i < h.attr.nodes; i++)
      {
        for (uint32_t j = 0; j <= i; j++)
        {
          if (j == i)
          {
            (*out) << param.fill << '\n';
          }
          else
          {
            if (e.index.i == i && e.index.j == j)
            {
              (*out) << e.scores[si].s << '\t';
              if (ctr < h.attr.edges)
              {
                e.unserialize(sf, h);
              }
            }
            else
            {
              (*out) << param.fill << '\t';
            }
          }
        }
      }
    }
    else // no diagonal
    {
      SeidrFileEdge e;
      e.unserialize(sf, h);
      uint64_t ctr = 1;
      for (uint32_t i = 1; i < h.attr.nodes; i++)
      {
        for (uint32_t j = 0; j < i; j++)
        {
          if (e.index.i == i && e.index.j == j)
          {
            (*out) << e.scores[si].s;
            if (ctr < h.attr.edges)
            {
              e.unserialize(sf, h);
            }
          }
          else
          {
            (*out) << param.fill;
          }
          *out << (j == i - 1 ? '\n' : '\t');
        }
      }
    }
  }
}

int adjacency(int argc, char * argv[]) {

  LOG_INIT_CERR();

  // Variables used by the function
  seidr_param_adjacency_t param;

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

    TCLAP::ValueArg<std::string>
    arg_outfmt("F", "outfmt", "Output format", false,
               "m", "m");
    cmd.add(arg_outfmt);

    TCLAP::ValueArg<std::string>
    arg_fill("m", "missing", "Filler character for missing nodes", false,
             "0", "0");
    cmd.add(arg_fill);

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

    TCLAP::SwitchArg
    switch_diag("D", "diagonal", "Print diagonal in triangular output", cmd,
                false);

    TCLAP::UnlabeledValueArg<std::string>
    arg_infile("infile", "Input file (aggregated gene counts)", true,
               "", "");

    cmd.add(arg_infile);

    // Parse arguments
    cmd.parse(argc, args);
    param.el_file = arg_infile.getValue();
    param.out_file = arg_outfile.getValue();
    param.tpos = arg_tindex.getValue();
    param.prec = arg_prec.getValue();
    param.force = switch_force.getValue();
    param.outfmt = arg_outfmt.getValue();
    param.fill = arg_fill.getValue();
    param.diag = switch_diag.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    std::cerr << "[Parser error]: " << except.what() << std::endl;
    return EINVAL;
  }

  try
  {
    file_can_read(param.el_file.c_str());

    if (param.out_file != "-")
    {
      if (file_exists(param.out_file) && ! param.force)
        throw std::runtime_error("File exists: " + param.out_file);

      file_can_create(param.out_file.c_str());
    }
  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << "[Runtime error]: " << except.what() << '\n';
    return errno;
  }

  if (param.outfmt == "m")
  {
    make_full(param);
  }
  else if (param.outfmt == "lm")
  {
    make_lower(param);
  }

  return 0;
}
