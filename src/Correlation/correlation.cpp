#include <cor_fun.h>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <common.h>
#include <fs.h>
#include <BSlogger.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>

int main(int argc, char ** argv)
{
  std::string infile;
  std::string outfile;
  std::string method;
  std::string gene_file;
  std::string targets_file;
  bool do_scale;
  bool abs;
  bool force = false;
  LOG_INIT_CLOG();
  unsigned verbosity;

  try
  {
    TCLAP::CmdLine cmd("Pearson/Spearman correlation for Seidr", ' ', version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<std::string>
    genefile_arg("g", "genes", "File containing gene names", true, "",
                 "string");
    cmd.add(genefile_arg);

    TCLAP::ValueArg<std::string>
    targets_arg("t", "targets", "File containing gene names"
                " of genes of interest. The network will only be"
                " calculated using these as the sources of potential connections.",
                false, "", "string");
    cmd.add(targets_arg);

    TCLAP::ValueArg<std::string>
    outfile_arg("o", "outfile", "Output file path", false, "edgelist.tsv",
                "edgelist.tsv");
    cmd.add(outfile_arg);

    std::vector<std::string> method_vec{"pearson", "spearman"};
    TCLAP::ValuesConstraint<std::string> allowed_methods( method_vec );
    TCLAP::ValueArg<std::string>
    method_arg("m", "method", "Correlation method <pearson>", false, "pearson",
               &allowed_methods);
    cmd.add(method_arg);

    TCLAP::ValueArg<unsigned>
    verbosity_arg("v", "verbosity", "Verbosity level (lower is less verbose)",
                  false, 3, "3");
    cmd.add(verbosity_arg);

    TCLAP::SwitchArg
    switch_scale("s", "scale", "Transform data to z-scores", cmd, false);

    TCLAP::SwitchArg
    switch_abs("a", "absolute", "Report absolute values", cmd, false);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.parse(argc, argv);
    infile = infile_arg.getValue();
    outfile = outfile_arg.getValue();
    do_scale = switch_scale.getValue();
    method = method_arg.getValue();
    abs = switch_abs.getValue();
    force = switch_force.getValue();
    gene_file = genefile_arg.getValue();
    targets_file = targets_arg.getValue();
    verbosity = verbosity_arg.getValue();

    log.set_log_level(verbosity);
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    log(LOG_ERR) << e.error() << " for arg " << e.argId() << '\n';
    return 1;
  }

  try
  {
    outfile = to_absolute(outfile);
    infile = to_absolute(infile);
    gene_file = to_absolute(gene_file);
    if (targets_file != "")
      targets_file = to_absolute(targets_file);

    if (! file_exists(dirname(outfile)) )
      throw std::runtime_error("Directory does not exist: " + dirname(outfile));

    if (! file_exists(infile) )
      throw std::runtime_error("File does not exist: " + infile);

    if (! file_exists(gene_file) )
      throw std::runtime_error("File does not exist: " + gene_file);

    if (! file_exists(targets_file) && targets_file != "" )
      throw std::runtime_error("File does not exist: " + targets_file);

    if (! regular_file(infile) )
      throw std::runtime_error("Not a regular file: " + infile);

    if (! file_can_read(infile) )
      throw std::runtime_error("Cannot read: " + infile);

    if (! force && file_exists(outfile))
      throw std::runtime_error("File exists: " + outfile);

  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << "[Runtime error]: " << e.what() << '\n';
    return 1;
  }

  arma::mat gm;

  try
  {
    gm.load(infile.c_str(), arma::raw_ascii);
    verify_matrix(gm);
    if (do_scale)
    {
      log(LOG_INFO) << "Transforming matrix to z-score\n";
      scale(gm);
    }

    if (method == "spearman")
    {
      log(LOG_INFO) << "Transforming matrix to value ranks\n";
      to_rank(gm);
    }

    gm = arma::cor(gm);
    if (targets_file == "")
    {
      write_lm(gm, outfile, abs);
    }
    else
    {
      std::vector<std::string> genes = read_genes(gene_file);
      std::vector<std::string> targets = read_genes(targets_file);
      std::unordered_map<std::string, uint64_t> gene_map;
      uint64_t ctr = 0;
      for (auto& g : genes)
        gene_map[g] = ctr++;
      std::ofstream ofs(outfile, std::ios::out);
      for (auto& t : targets)
      {
        uint64_t i;
        try
        {
          i = gene_map.at(t);
        }
        catch (std::exception& e)
        {
          log(LOG_ERR) << e.what() << '\n';
          log(LOG_ERR) << "Target gene " << t << "is not in the expression matrix\n";
        }
        for (uint64_t j = 0; j < gm.n_cols; j++)
        {
          if (i == j) continue;
          ofs << genes[i] << '\t' << genes[j] << '\t' <<
              (abs ? fabs(gm(i, j)) : gm(i, j)) << '\n';
        }
      }
    }
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << "[Runtime error]: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
