// Seidr
#include <common.h>
#include <compare.h>
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
#include <boost/numeric/conversion/cast.hpp>

cindex::cindex(SeidrFileHeader& a, SeidrFileHeader& b)
{
  std::set<std::string> set_a;
  uint32_t ctr = 0;
  for (auto& gene : a.nodes)
  {
    map_a[gene] = ctr++;
    set_a.insert(gene);
    uni.insert(gene);
  }
  // reset counter
  ctr = 0;
  for (auto& gene : b.nodes)
  {
    map_b[gene] = ctr++;
    if (set_a.find(gene) != set_a.end())
    {
      intersect.insert(gene);
    }
    else
    {
      diff_b.insert(gene);
    }
    uni.insert(gene);
  }
  for (auto& gene : a.nodes)
  {
    if (intersect.find(gene) == intersect.end())
    {
      diff_a.insert(gene);
    }
  }
  ctr = 0;
  for (auto& gene : uni)
  {
    map_c[gene] = ctr++;
  }
  // Create vector of unique merged node names
  for (auto& n : uni)
    node_uni.push_back(n);
}

cindex::cindex(SeidrFileHeader& a, SeidrFileHeader& b, std::string& file_orth)
{
  orth_mode = true;
  std::ifstream ifs_orth(file_orth, std::ios::in);
  std::string line;
  while (std::getline(ifs_orth, line))
  {
    std::string left, right;
    std::stringstream ss(line);
    ss >> left;
    ss >> right;
    dict_orth[left] = right;
    dict_orth[right] = left;
  }

  std::set<std::string> set_a;
  uint32_t ctr = 0;
  for (auto& gene : a.nodes)
  {
    map_a[gene] = ctr++;
    set_a.insert(gene);
    uni.insert(gene);
  }
  // reset counter
  std::set<std::string> set_b;
  ctr = 0;
  for (auto& gene : b.nodes)
  {
    map_b[gene] = ctr++;
    set_b.insert(gene);
    uni.insert(gene);
  }

  ctr = 0;
  for (auto& gene : uni)
  {
    auto ptr = dict_orth.find(gene);
    if (ptr != dict_orth.end())
    {
      if (map_c.find(ptr->first) == map_c.end())
      {
        map_c[ptr->first] = ctr;
        map_c[ptr->second] = ctr;
        dict_index[ctr] = string_p(ptr->first, ptr->second);
        ctr++;
      }
    }
    else
    {
      map_c[gene] = ctr;
      ctr++;
    }
  }
  // Create vector of unique merged node names
  node_uni.resize(ctr);
  for (auto& n : map_c)
  {
    if (dict_index.find(n.second) != dict_index.end())
    {
      std::string s1 = dict_index[n.second].first;
      std::string s2 = dict_index[n.second].second;
      node_uni[n.second] = s1 + "," + s2;
    }
    else
    {
      node_uni[n.second] = n.first;
    }
  }
}


uint32_p remap_index(SeidrFileEdge& e, SeidrFileHeader& h, cindex& index)
{
  // Need to copy from SeidrFileEdge& e since we can't refer to fields of a
  // packed struct.
  uint32_t i = e.index.i;
  uint32_t j = e.index.j;
  uint32_p p1(i, j);
  uint32_p p1r = index.reindex(h, e);
  return p1r;
}

SeidrFileEdge create_edge(SeidrFileEdge& e1, SeidrFileHeader& h1,
                          SeidrFileEdge& e2, SeidrFileHeader& h2, uint32_p& p,
                          uint8_t score, uint32_t& ti1, uint32_t& ti2)
{
  SeidrFileEdge e;
  e.index.i = p.first;
  e.index.j = p.second;
  e.supp_int.push_back(score);
  e.supp_flt = {0, 0};
  edge_score es;
  if (score == 0)
  {
    es.r = e1.scores[ti1].r + e2.scores[ti2].r;
    e.supp_flt[0] = boost::numeric_cast<float>(e1.scores[ti1].r);
    e.supp_flt[1] = boost::numeric_cast<float>(e2.scores[ti2].r);
  }
  else if (score == 1)
  {
    es.r = e1.scores[ti1].r;
    e.supp_flt[0] = boost::numeric_cast<float>(e1.scores[ti1].r);
  }
  else
  {
    es.r = e2.scores[ti2].r;
    e.supp_flt[1] = boost::numeric_cast<float>(e2.scores[ti2].r);
  }
  EDGE_SET_EXISTING(e.attr.flag);
  e.scores.push_back(es);
  return e;
}

void finalize(std::string& file_out, SeidrFileHeader& h,
              double min, double max)
{
  std::string tmp = file_out + ".tmp";
  SeidrFile out(tmp.c_str());
  out.open("r");
  SeidrFileHeader ho;
  ho.unserialize(out);

  SeidrFile fin(file_out.c_str());
  fin.open("w");
  h.serialize(fin);

  if (h.attr.dense)
  {

  }
  else
  {
    for (uint64_t i = 0; i < h.attr.edges; i++)
    {
      SeidrFileEdge e;
      e.unserialize(out, ho);
      e.scores[0].s = unity_stand(min, max, e.scores[0].r);
      e.serialize(fin, h);
    }
  }

  out.close();
  fin.close();

  remove(tmp);

}

void compare_nodes(SeidrFile& f1, SeidrFileHeader& h1,
                   SeidrFile& f2, SeidrFileHeader& h2,
                   std::string& trans)
{
  cindex index;
  if (trans == "")
    index = cindex(h1, h2);
  else
    index = cindex(h1, h2, trans);

  std::vector<std::string> node_names = index.node_union();
  std::map<std::string, int> nm;
  for (std::string& n : node_names)
    nm[n] = 3;
  f1.seek(0);
  f1.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    uint32_p re = index.reindex(h, e);
    nm[node_names[re.first]] = 1;
    nm[node_names[re.second]] = 1;
  });
  f2.seek(0);
  f2.each_edge([&](SeidrFileEdge & e, SeidrFileHeader & h) {
    uint32_p re = index.reindex(h, e);
    auto ptr = nm.find(node_names[re.first]);
    if (ptr->second == 1 || ptr->second == 0)
      ptr->second = 0;
    else
      ptr->second = 2;

    ptr = nm.find(node_names[re.second]);
    if (ptr->second == 1 || ptr->second == 0)
      ptr->second = 0;
    else
      ptr->second = 2;
  });
  for (auto& n : nm)
  {
    std::cout << n.first << '\t' << n.second << '\n';
  }
  f1.close();
  f2.close();
}

void next_edge(SeidrFileEdge& e, SeidrFile& f, SeidrFileHeader& h,
               uint64_t& i, uint64_t& j)
{
  if (h.attr.dense)
  {
    e.unserialize(f, h);
    e.index.i = i; e.index.j = j;
    if (j == (i - 1))
    {
      i++;
      j = 0;
    }
    else
    {
      j++;
    }
    if (! EDGE_EXISTS(e.attr.flag))
      next_edge(e, f, h, i, j);
  }
  else
  {
    e.unserialize(f, h);
  }
}

void compare_ss(SeidrFile& f1, SeidrFileHeader& h1,
                SeidrFile& f2, SeidrFileHeader& h2,
                std::string& file_out,
                uint32_t tindex1, uint32_t tindex2,
                SeidrFileHeader& h3, std::string& trans)
{
  std::string tmp = file_out + ".tmp";
  SeidrFile out(tmp.c_str());
  out.open("w");

  uint64_t e1i = 1;
  uint64_t e1j = 0;
  uint64_t e2i = 1;
  uint64_t e2j = 0;

  cindex index;
  if (trans == "")
    index = cindex(h1, h2);
  else
    index = cindex(h1, h2, trans);

  SeidrFileEdge e1;
  next_edge(e1, f1, h1, e1i, e1j);

  SeidrFileEdge e2;
  next_edge(e2, f2, h2, e2i, e2j);

  uint32_p p1r = remap_index(e1, h1, index);
  uint32_p p2r = remap_index(e2, h2, index);

  uint64_t ctr1 = 1;
  uint64_t ctr2 = 1;
  uint64_t ctr3 = 0;

  bool break1 = false;
  bool break2 = false;

  std::vector<std::string>& node_names = index.node_union();

  h3.attr.nodes = node_names.size();
  h3.attr.edges = 8;
  h3.attr.nalgs = 1;
  h3.attr.dense = 0;
  h3.attr.nsupp = 3;
  h3.attr.nsupp_flt = 2;
  h3.attr.nsupp_int = 1;
  h3.supp = {"N", "A", "B"};
  h3.algs = {"Merge"};

  h3.nodes = node_names;

  h3.serialize(out);

  double min = std::numeric_limits<double>::infinity();
  double max = -std::numeric_limits<double>::infinity();

  while (true)
  {
    SeidrFileEdge e3;
    if (p1r.first < p2r.first)
    {
      e3 = create_edge(e1, h1, e2, h2, p1r, 1, tindex1, tindex2);
      e3.serialize(out, h3);
      // Check if this was the last edge in the file
      if (break1)
        break;
      next_edge(e1, f1, h1, e1i, e1j);
      p1r = remap_index(e1, h1, index);
      ctr1++;
      if (ctr1 == h1.attr.edges)
        break1 = true;
    }
    else if (p1r.first > p2r.first)
    {
      e3 = create_edge(e1, h1, e2, h2, p2r, 2, tindex1, tindex2);
      e3.serialize(out, h3);
      // Check if this was the last edge in the file
      if (break2)
        break;
      next_edge(e2, f2, h2, e2i, e2j);
      p2r = remap_index(e2, h2, index);
      ctr2++;
      if (ctr2 == h2.attr.edges)
        break2 = true;
    }
    else
    {
      if (p1r.second < p2r.second)
      {
        e3 = create_edge(e1, h1, e2, h2, p1r, 1, tindex1, tindex2);
        e3.serialize(out, h3);
        // Check if this was the last edge in the file
        if (break1)
          break;
        next_edge(e1, f1, h1, e1i, e1j);
        p1r = remap_index(e1, h1, index);
        ctr1++;
        if (ctr1 == h1.attr.edges)
          break1 = true;
      }
      else if (p1r.second > p2r.second)
      {
        e3 = create_edge(e1, h1, e2, h2, p2r, 2, tindex1, tindex2);
        e3.serialize(out, h3);
        if (break2)
          break;
        next_edge(e2, f2, h2, e2i, e2j);
        p2r = remap_index(e2, h2, index);
        ctr2++;
        if (ctr2 == h2.attr.edges)
          break2 = true;
      }
      else
      {
        e3 = create_edge(e1, h1, e2, h2, p1r, 0, tindex1, tindex2);
        e3.serialize(out, h3);
        if (break1 && break2)
        {
          break;
        }
        else if (break1)
        {
          next_edge(e2, f2, h2, e2i, e2j);
          p2r = remap_index(e2, h2, index);
          ctr2++;
          ctr3++;
          break;
        }
        else if (break2)
        {
          next_edge(e1, f1, h1, e1i, e1j);
          p1r = remap_index(e1, h1, index);
          ctr1++;
          ctr3++;
          break;
        }
        else
        {
          next_edge(e1, f1, h1, e1i, e1j);
          p1r = remap_index(e1, h1, index);
          next_edge(e2, f2, h2, e2i, e2j);
          p2r = remap_index(e2, h2, index);
        }

        ctr1++;
        ctr2++;
        ctr3++;
        if (ctr1 == h1.attr.edges)
          break1 = true;
        if (ctr2 == h2.attr.edges)
          break2 = true;
      }
    }
    min = e3.scores[0].r < min ? e3.scores[0].r : min;
    max = e3.scores[0].r > max ? e3.scores[0].r : max;
  }
  SeidrFileEdge e3;
  if (break1 && !break2)
  {
    for (uint64_t i = ctr2 - 1; i < h2.attr.edges; i++)
    {
      e3 = create_edge(e1, h1, e2, h2, p2r, 2, tindex1, tindex2);
      e3.serialize(out, h3);
      if (i < h2.attr.edges - 1)
      {
        next_edge(e2, f2, h2, e2i, e2j);
        p2r = remap_index(e2, h2, index);
      }
      min = e3.scores[0].r < min ? e3.scores[0].r : min;
      max = e3.scores[0].r > max ? e3.scores[0].r : max;
      ctr2++;
    }
  }
  else if (break2 && !break1)
  {
    for (uint64_t i = ctr1 - 1; i < h1.attr.edges; i++)
    {
      e3 = create_edge(e1, h1, e2, h2, p1r, 1, tindex1, tindex2);
      e3.serialize(out, h3);
      if (i < h1.attr.edges - 1)
      {
        next_edge(e1, f1, h1, e1i, e1j);
        p1r = remap_index(e1, h1, index);
      }
      min = e3.scores[0].r < min ? e3.scores[0].r : min;
      max = e3.scores[0].r > max ? e3.scores[0].r : max;
      ctr1++;
    }
  }
  out.close();
  // Calc standardized scores
  uint64_t ne = ctr1 + ctr2 - ctr3 - 1;
  h3.attr.edges = ne;
  finalize(file_out, h3, min, max);
}

int compare(int argc, char * argv[]) {

  logger log(std::cerr, "compare");

  // Variables used by the function
  std::string net1;
  std::string net2;
  std::string file_out;
  uint32_t tindex1;
  uint32_t tindex2;
  std::string trans;
  bool node_mode;
  bool force = false;

  // We ignore the first argument
  const char * args[argc - 1];
  std::string pr(argv[0]);
  pr += " compare";
  args[0] = pr.c_str();
  for (int i = 2; i < argc; i++) args[i - 1] = argv[i];
  argc--;

  try
  {
    // Add arguments from the command line
    TCLAP::CmdLine cmd("Compare edges or nodes in two networks.\n", 
                       ' ', version);

    TCLAP::ValueArg<std::string>
    arg_out("o", "outfile", "Output file (Merged)", false, "compare.sf", 
            "string");
    cmd.add(arg_out);

    TCLAP::ValueArg<uint32_t>
    arg_tindex1("i", "index-a", "Merge scores on this index for network A", 
                false, 0, "int");

    cmd.add(arg_tindex1);
    TCLAP::ValueArg<uint32_t>
    arg_tindex2("j", "index-b", "Merge scores on this index for network B", 
                false, 0, "int");

    cmd.add(arg_tindex2);

    TCLAP::ValueArg<std::string>
    arg_trans("t", "translate", "Translate node names in network 1 according "
              "to this table", false, "", "string");
    cmd.add(arg_trans);

    TCLAP::SwitchArg
    switch_nodes("n", "nodes", "Print node overlap", cmd, false);

    TCLAP::UnlabeledValueArg<std::string>
    arg_net1("network-1",
             "Input file (aggregated gene counts)", true, "", "string");
    cmd.add(arg_net1);

    TCLAP::UnlabeledValueArg<std::string>
    arg_net2("network-2",
             "Input file (aggregated gene counts)", true, "", "string");

    TCLAP::SwitchArg
    switch_force("f", "force", "Force overwrite if output already exists", cmd,
                 false);

    cmd.add(arg_net2);

    // Parse arguments
    cmd.parse(argc, args);
    net1 = arg_net1.getValue();
    net2 = arg_net2.getValue();
    file_out = arg_out.getValue();
    tindex1 = arg_tindex1.getValue();
    tindex2 = arg_tindex2.getValue();
    node_mode = switch_nodes.getValue();
    trans = arg_trans.getValue();
    force = switch_force.getValue();
  }
  catch (TCLAP::ArgException& except)
  {
    log(LOG_ERR) << except.what() << '\n';
    return 1;
  }

  try
  {
    net1 = to_canonical(net1);
    net2 = to_canonical(net2);
    if (file_exists(file_out) && ! force)
      throw std::runtime_error("File exists: " + file_out);
  }
  catch (std::runtime_error& e)
  {
    log(LOG_ERR) << e.what() << '\n';
    return errno;
  }

  SeidrFile n1(net1.c_str());
  n1.open("r");

  SeidrFile n2(net2.c_str());
  n2.open("r");

  SeidrFileHeader h1;
  h1.unserialize(n1);

  SeidrFileHeader h2;
  h2.unserialize(n2);

  if (tindex1 == 0)
    tindex1 = h1.attr.nalgs - 1;
  else
    tindex1--;

  if (tindex2 == 0)
    tindex2 = h2.attr.nalgs - 1;
  else
    tindex2--;

  if (node_mode)
  {
    compare_nodes(n1, h1, n2, h2, trans);
  }
  else
  {
    SeidrFileHeader h3;
    h3.version_from_char(_XSTR(VERSION));
    h3.cmd_from_args(argv, argc);
    compare_ss(n1, h1, n2, h2, file_out, tindex1, tindex2, h3, trans);
  }

  n1.close();
  n2.close();

  return 0;
}
