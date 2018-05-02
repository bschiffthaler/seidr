#include <pcor-fun.h>
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
  bool do_scale = true;
  bool abs;
  bool force = false;
  LOG_INIT_CLOG();

  try
  {
    TCLAP::CmdLine cmd("Partial correlation for Seidr", ' ', version);

    TCLAP::ValueArg<std::string>
    infile_arg("i", "infile", "The expression table (without headers)", true,
               "", "");
    cmd.add(infile_arg);

    TCLAP::ValueArg<std::string>
    outfile_arg("o", "outfile", "Output file path", false, "edgelist.tsv",
                "edgelist.tsv");
    cmd.add(outfile_arg);

    TCLAP::SwitchArg
    switch_abs("a", "absolute", "Report absolute values", cmd, false);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.parse(argc, argv);
    infile = infile_arg.getValue();
    outfile = outfile_arg.getValue();
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

    if(! force && file_exists(outfile))
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
    arma::mat pc = pcor(gm);
    write_lm(pc, outfile, abs);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << "[Runtime error]: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
