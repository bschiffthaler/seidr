#include <cor_fun.h>
#include <armadillo>
#include <tclap/CmdLine.h>
#include <common.h>
#include <fs.h>
#include <BSlogger.h>
#include <vector>
#include <string>

int main(int argc, char ** argv)
{
  std::string infile;
  std::string outfile;
  std::string method;
  bool do_scale;
  bool abs;
  bool force = false;
  LOG_INIT_CLOG();

  try
  {
    TCLAP::CmdLine cmd("Pearson/Spearman correlation for Seidr", ' ', version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

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

    if (! file_exists(dirname(outfile)) )
      throw std::runtime_error("Directory does not exist: " + dirname(outfile));

    if (! file_exists(infile) )
      throw std::runtime_error("File does not exist: " + infile);

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
    write_lm(gm, outfile, abs);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << "[Runtime error]: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
