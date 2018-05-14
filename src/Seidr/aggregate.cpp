/**
 * @file
 * @author Bastian Schiffthaler <bastian.schiffthaler@umu.se>
 *
 * @param infs Input files to be aggregated in binary format from el2bin
 * @param genes_file File containing names of genes present in the network
 * @param method Aggregation method: 'borda' implements borda counting
 *               'last' ranks missing edges last, 'ignore' ignores missing.
 * @return int 0 if the function succeeded, an error code otherwise.
 *
 * @section DESCRIPTION
 *
 * This function aggregates any number of networks (represented as binary
 * files by the el2bin function) into a meta-network of scores in [0, 1].
 * If the chosen method is borda counting, we calculate the mean rank of
 * each edge, ranking missing vertices last. The top1 algorithm presents
 * the final score of an edge as the best score of any algorithm. The
 * neural algorithm learns how to best predict new edges from a set of
 * gold standard gene interactions.
 */

// Seir
#include <common.h>
#include <Serialize.h>
#include <fs.h>
#include <BSlogger.h>
// External
#include <vector>
#include <cmath>
#include <functional>
#include <string>
#include <fstream>
#include <forward_list>
#include <algorithm>
#include <tclap/CmdLine.h>
#include <stdexcept>
#include <cerrno>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>

using boost::lexical_cast;
using boost::numeric_cast;
namespace fs = boost::filesystem;

struct aggr_rank_t {
  double score = 0;
  double rank = 0;
  uint64_t index = 0;
};

bool ev_score_sort(aggr_rank_t a, aggr_rank_t b) { return (a.score < b.score);}
bool ev_index_sort(aggr_rank_t a, aggr_rank_t b) { return (a.index < b.index);}

void rank_vector(std::vector<aggr_rank_t>& ev)
{
  logger log(std::cerr, "rank_vector");
  std::sort(ev.begin(), ev.end(), ev_score_sort);
  log(LOG_DEBUG) << "Computing ranks\n";
  auto it = ev.begin();
  uint64_t pos = 0;
  double prev = it->score;
  uint64_t start = 0;
  double rank;
  while (it != ev.end())
  {
    it++; pos++;
    if (! almost_equal(it->score, prev) || it == ev.end())
    {
      rank = ( lexical_cast<double>(pos) + 1 +
               lexical_cast<double>(start) ) / 2;
      for (size_t i = start; i < pos; i++)
      {
        ev[i].rank = rank;
      }
      if (it != ev.end())
      {
        start = pos;
        prev = it->score;
      }
    }
  }
  std::sort(ev.begin(), ev.end(), ev_index_sort);
}

void aggr_top1(std::vector<SeidrFileHeader>& header_vec,
               std::vector<SeidrFileEdge>& v,
               SeidrFileEdge& result,
               std::vector<uint8_t>& existing,
               double& rank,
               uint64_t& ex_sum)
{
  uint16_t index = 0;
  for (uint16_t index_b = 0; index_b < result.scores.size(); index_b++)
  {
    if (existing[index_b] && result.scores[index_b].r < rank)
    {
      rank = result.scores[index_b].r;
      index = index_b;
    }
  }

  edge_score es;
  es.r = rank;
  es.s = 0;
  result.scores.push_back(es);

  result.supp_int.push_back(numeric_cast<int>(index));

  if (ex_sum > 0)
  {
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_DIRECT(v[index].attr.flag))
    {
      if (EDGE_IS_AB(v[index].attr.flag))
      {
        EDGE_SET_AB(result.attr.flag);
      }
      else
      {
        EDGE_SET_BA(result.attr.flag);
      }
    }
  }
}

void aggr_top2(std::vector<SeidrFileHeader>& header_vec,
               std::vector<SeidrFileEdge>& v,
               SeidrFileEdge& result,
               std::vector<uint8_t>& existing,
               double& rank,
               uint64_t& ex_sum)
{

  uint16_t index1 = 0;
  uint16_t index2 = 0;
  double r1 = rank;
  double r2 = rank;
  for (uint16_t index_b = 0; index_b < result.scores.size(); index_b++)
  {
    if (existing[index_b] && result.scores[index_b].r < r1)
    {
      r1 = result.scores[index_b].r;
      index1 = index_b;
    }
    else if (existing[index_b] && result.scores[index_b].r < r2)
    {
      r2 = result.scores[index_b].r;
      index2 = index_b;
    }
  }

  edge_score es;
  es.r = (r1 + r2) / 2;
  es.s = 0;
  result.scores.push_back(es);

  result.supp_int.push_back(numeric_cast<int>(index1));
  result.supp_int.push_back(numeric_cast<int>(index2));

  if (ex_sum > 0)
  {
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_DIRECT(v[index1].attr.flag) && EDGE_IS_DIRECT(v[index2].attr.flag))
    {
      if (EDGE_IS_AB(v[index1].attr.flag) && EDGE_IS_AB(v[index2].attr.flag))
      {
        EDGE_SET_AB(result.attr.flag);
      }
      else if (EDGE_IS_BA(v[index1].attr.flag) && EDGE_IS_BA(v[index2].attr.flag))
      {
        EDGE_SET_BA(result.attr.flag);
      }
    }
  }
}

void aggr_borda(std::vector<SeidrFileHeader>& header_vec,
                std::vector<SeidrFileEdge>& v,
                SeidrFileEdge& result,
                std::vector<uint8_t>& existing,
                double& rank,
                uint64_t& ex_sum)
{
  uint16_t index = 0;
  double sum = 0;
  double cnt = 0;
  for (uint16_t index_b = 0; index_b < result.scores.size(); index_b++)
  {
    if (existing[index_b])
    {
      sum += result.scores[index_b].r;
      cnt++;
    }
  }

  edge_score es;
  es.r = sum / cnt;
  es.s = 0;
  result.scores.push_back(es);

  if (ex_sum > 0)
  {
    EDGE_SET_EXISTING(result.attr.flag);
    if (EDGE_IS_DIRECT(v[index].attr.flag))
    {
      if (EDGE_IS_AB(v[index].attr.flag))
      {
        EDGE_SET_AB(result.attr.flag);
      }
      else
      {
        EDGE_SET_BA(result.attr.flag);
      }
    }
  }
}

void aggr_irp(std::vector<SeidrFileHeader>& header_vec,
              std::vector<SeidrFileEdge>& v,
              SeidrFileEdge& result,
              std::vector<uint8_t>& existing,
              double& rank,
              uint64_t& ex_sum)
{
  uint16_t index = 0;
  double r;
  double nodes = header_vec[0].attr.nodes;
  double nmax = log10((nodes * (nodes - 1)) / 2);
  uint16_t cnt = 0;
  if (existing[0])
  {
    // As LHS and RHS are logs multiplication of both
    // becomes a sum
    r = 1 + log10(result.scores[0].r);
    cnt++;
  }
  else
  {
    r = 1 + nmax;
  }
  for (uint16_t index_b = 1; index_b < result.scores.size(); index_b++)
  {
    if (existing[index_b])
    {
      r += (1 + log10(result.scores[index_b].r));
      cnt++;
    }
    else
    {
      r += (1 + nmax);
    }
  }

  edge_score es;
  es.r = r;
  es.s = 0;
  result.scores.push_back(es);

  if (cnt > 0)
  {
    EDGE_SET_EXISTING(result.attr.flag);
  }
}


SeidrFileEdge calc_score(std::vector<SeidrFileHeader>& header_vec,
                         std::vector<SeidrFileEdge>& v,
                         std::vector<uint8_t>& get_next,
                         uint32_t i, uint64_t j,
                         std::function<void(std::vector<SeidrFileHeader>&,
                             std::vector<SeidrFileEdge>&,
                             SeidrFileEdge&,
                             std::vector<uint8_t>&,
                             double&,
                             uint64_t&)> aggr_fun)
{
  SeidrFileEdge result;
  result.index.i = i;
  result.index.j = j;
  std::vector<uint8_t> existing;
  uint64_t ex_sum = 0;
  double rank = std::numeric_limits<double>::infinity();
  uint32_t index_a = 0;
  edge_score dummy;
  dummy.r = std::numeric_limits<double>::quiet_NaN();
  dummy.s = std::numeric_limits<double>::quiet_NaN();
  for (auto& e : v)
  {
    if (e.index.i == i && e.index.j == j)
    {
      if (EDGE_EXISTS(e.attr.flag))
      {
        for (auto& s : e.scores)
        {
          existing.push_back(1);
          ex_sum += 1;
          result.scores.push_back(s);
        }
        for (auto& s : e.supp_str)
          result.supp_str.push_back(s);
        for (auto& s : e.supp_int)
          result.supp_int.push_back(s);
        for (auto& s : e.supp_flt)
          result.supp_flt.push_back(s);
      }
      else
      {
        existing.push_back(0);
        result.scores.push_back(dummy);
        for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_str; a++)
          result.supp_str.push_back("");
        for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_int; a++)
          result.supp_int.push_back(0);
        for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_flt; a++)
          result.supp_flt.push_back(0);
      }
      get_next[index_a] = 1;
    }
    else
      // Storage is sparse and rank file is ahead of current edge
      // We fill the edge with dummy data
    {
      for (uint16_t a = 0; a < header_vec[index_a].attr.nalgs; a++)
      {
        result.scores.push_back(dummy);
        existing.push_back(0);
      }
      for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_str; a++)
        result.supp_str.push_back("");
      for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_int; a++)
        result.supp_int.push_back(0);
      for (uint16_t a = 0; a < header_vec[index_a].attr.nsupp_flt; a++)
        result.supp_flt.push_back(0);
      // Do not read next edge from this file
      get_next[index_a] = 0;
    }
    index_a++;
  }

  // Call function pointer
  aggr_fun(header_vec, v, result, existing, rank, ex_sum);

  return result;
}

int aggregate(int argc, char * argv[])
{

  logger log(std::cerr, "aggregate");

  // Variables used by the function
  std::vector<std::string> infs;
  std::string method;
  std::string out_file;
  bool force = false;

  // We ignore the first argument, the function name
#ifndef TEST_BUILD
  const char * args[argc - 1];
  std::string a0(argv[0]);
  a0 += " aggregate";
  const char * a1 = a0.c_str();
  args[0] = a1;
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;
#else
  char ** args = argv;
#endif

  std::vector<std::string> method_constraints{"top1", "top2", "borda", "irp"};

  try
  {
    TCLAP::ValuesConstraint<std::string> constraints_m(method_constraints);
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Aggregate a set of ranked representation networks.", 
                       ' ', version);

    TCLAP::ValueArg<std::string> 
    arg_method("m", "method", "Method to aggregate networks <irp>.", false,
                                            "irp", &constraints_m);
    cmd.add(arg_method);

    TCLAP::ValueArg<std::string> 
    arg_outfile("o", "outfile", "Where to write the output file.", false,
        "aggregated.sf", "aggregated.sf");
    cmd.add(arg_outfile);

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    TCLAP::UnlabeledMultiArg<std::string> 
    arg_input_files("infiles", "Input files", true, "", "");

    cmd.add(arg_input_files);

    // Parse arguments
    cmd.parse(argc, args);
    method = arg_method.getValue();
    infs = arg_input_files.getValue();
    out_file = std::string(arg_outfile.getValue());
    force = switch_force.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << "[Invalid argument exception]: " << except.what() << '\n';
    return EINVAL;
  }

  try
  {
    log(LOG_INFO) << out_file << '\n';
    out_file = to_absolute(out_file);

    if (file_exists(out_file) && ! force)
      throw std::runtime_error("File exists: " + out_file);

    if (! file_can_create(out_file))
      throw std::runtime_error("Can't create file " + out_file + "\n");

    for (size_t i = 0; i < infs.size(); i++)
    {
      infs[i] = to_canonical(infs[i]);
      if (! file_can_read(infs[i]))
        throw std::runtime_error("Can't read file " + infs[i] + "\n");
      log(LOG_INFO) << "Have file " << infs[i] << '\n';
    }

  }
  catch (std::runtime_error& except)
  {
    log(LOG_ERR) << "[Runtime exception]: " << except.what() << '\n';
    return errno;
  }

  std::vector<SeidrFile> infile_vec;
  std::vector<SeidrFileHeader> header_vec;

  log(LOG_INFO) << "Getting headers\n";
  uint32_t counter = 0;
  for (auto& f : infs)
  {
    // Create a SeidrFile from input string and call bgzf_open()
    infile_vec.push_back(SeidrFile(f.c_str()));
    infile_vec[counter].open("r");
    // Read file header and push into vector
    SeidrFileHeader h;
    h.unserialize(infile_vec[counter]);
    header_vec.push_back(h);
    counter++;
  }

  log(LOG_DEBUG) << "Create new header\n";
  // Setup a new header (temporary) for serializing
  SeidrFileHeader h;
  h.attr.nodes = header_vec[0].attr.nodes;
  h.attr.dense = 1;
  h.attr.nalgs = infile_vec.size();
  for (auto& hh : header_vec)
  {
    h.attr.nsupp += hh.attr.nsupp;
    h.attr.nsupp_str += hh.attr.nsupp_str;
    h.attr.nsupp_int += hh.attr.nsupp_int;
    h.attr.nsupp_flt += hh.attr.nsupp_flt;
    for (auto s : hh.supp)
    {
      h.supp.push_back(s);
    }
    for (auto a : hh.algs)
    {
      h.algs.push_back(a);
    }
  }
  h.nodes = header_vec[0].nodes;
  h.version_from_char(_XSTR(VERSION));
  h.cmd_from_args(argv, argc);
  h.algs.push_back(method);
  h.attr.nalgs += 1;

  if (method == "top1")
  {
    h.attr.nsupp += 1;
    h.attr.nsupp_int += 1;
    h.supp.push_back("T");
  }

  if (method == "top2")
  {
    h.attr.nsupp += 2;
    h.attr.nsupp_int += 2;
    h.supp.push_back("T1");
    h.supp.push_back("T2");
  }


  log(LOG_INFO) << "Aggregating using method: " << method << '\n';
  std::string tempfile = out_file + std::string(".tmp");
  tempfile = to_absolute(tempfile);
  log(LOG_INFO) << "Using temp file: " << tempfile << '\n';
  SeidrFile tmp(tempfile.c_str());
  tmp.open("w");
  h.serialize(tmp);
  double min = std::numeric_limits<double>::infinity();
  double max = -std::numeric_limits<double>::infinity();
  // Aggregate in lower triangular order (guaranteed from `seidr import`)
  // First, setup vectors of length equal to the number of input files
  // One holds the edges in the input files, the other decides if we need
  // to read the next edge from this file
  std::vector<SeidrFileEdge> v;
  std::vector<uint8_t> get_next;
  std::vector<uint64_t> edges_read;
  v.resize(counter);
  get_next.resize(counter);
  edges_read.resize(counter);
  for (auto& e : edges_read) e = 0;
  // Define aggregation method
  std::function<void(std::vector<SeidrFileHeader>&,
                     std::vector<SeidrFileEdge>&,
                     SeidrFileEdge&,
                     std::vector<uint8_t>&,
                     double&,
                     uint64_t&)> aggr_fun;
  if (method == "top1")
    aggr_fun = aggr_top1;
  else if (method == "borda")
    aggr_fun = aggr_borda;
  else if (method == "top2")
    aggr_fun = aggr_top2;
  else if (method == "irp")
    aggr_fun = aggr_irp;
  else
    throw std::runtime_error("Unknown ranking method: " + method);
  // In the first run, we read once from all files
  for (auto& x : get_next)
    x = 1;
  uint64_t index = 0;
  std::vector<aggr_rank_t> rvec;
  for (uint64_t i = 1; i < h.attr.nodes; i++)
  {
    for (uint64_t j = 0; j < i; j++)
    {
      for (uint32_t k = 0; k < counter; k++)
      {
        if (get_next[k])
        {
          edges_read[k]++;
          // If the last edge in a sparse network was read, do not attempt a
          // file read
          if (! header_vec[k].attr.dense &&
              edges_read[k] > header_vec[k].attr.edges)
          {
            get_next[k] = 0;
            SeidrFileEdge e;
            e.index.i = std::numeric_limits<uint32_t>::infinity();
            e.index.j = std::numeric_limits<uint32_t>::infinity();
            v[k] = e;
          }
          else
          {
            SeidrFileEdge e;
            e.unserialize(infile_vec[k], header_vec[k]);
            if (header_vec[k].attr.dense)
            {
              e.index.i = i;
              e.index.j = j;
            }
            v[k] = e;
          }
        }

      }
      SeidrFileEdge res = calc_score(header_vec, v, get_next, i, j, aggr_fun);
      if (EDGE_EXISTS(res.attr.flag))
      {
        h.attr.edges++;
        aggr_rank_t r;
        r.index = index;
        r.score = res.scores[counter].r;
        rvec.push_back(r);
        if (res.scores[counter].r < min)
          min = res.scores[counter].r;
        if (res.scores[counter].r > max)
          max = res.scores[counter].r;
      }
      res.serialize(tmp, h);
      index++;
      if (index % 100000 == 0)
        log(LOG_INFO) << "Processed " << index << " edges\n";
    }
  }
  log(LOG_INFO) << "Done aggregating. Ranking and standardizing aggregated edges\n";
  // call bgzf_close() on all open files
  for (auto& f : infile_vec)
    f.close();
  tmp.close();
  tmp.open("r");

  rank_vector(rvec);

  SeidrFileHeader ha;
  ha.unserialize(tmp);

  SeidrFile out(out_file.c_str());
  out.open("w");

  // Determine storage model of final file
  double ne = h.attr.edges;
  double nn = h.attr.nodes;
  if (ne < ((nn * (nn - 1) / 2) * 0.66))
    h.attr.dense = 0;
  // Write header
  h.serialize(out);

  // Serialize all edges
  index = 0;
  for (uint64_t i = 1; i < h.attr.nodes; i++)
  {
    for (uint64_t j = 0; j < i; j++)
    {
      SeidrFileEdge e;
      e.unserialize(tmp, ha);
      if (EDGE_EXISTS(e.attr.flag))
      {
        log(LOG_DEBUG) << rvec[index].index << ", "
                       << e.scores[counter].r << ", "
                       << rvec[index].score << ", "
                       << index << '\n';
        e.scores[counter].s = unity_stand(min, max, e.scores[counter].r);
        e.scores[counter].r = rvec[index].rank;
        index++;
      }
      if (h.attr.dense)
      {
        e.serialize(out, h);
      }
      else if (EDGE_EXISTS(e.attr.flag))
      {
        e.index.i = i;
        e.index.j = j;
        e.serialize(out, h);
      }
    }
  }

  tmp.close();
  out.close();

  log(LOG_INFO) << "Removing temp file: " << tempfile << '\n';
  fs::remove(fs::path(tempfile));
  log(LOG_INFO) << "Finished\n";
  return 0;
}

#ifdef TEST_BUILD
int main(int argc, char *argv[])
{
  aggregate(argc, argv);
  return 0;
}
#endif
