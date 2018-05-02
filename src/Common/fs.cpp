#include <fstream>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

bool exists(const std::string fname)
{
  fs::path p(fname);
  if(! fs::exists(p))
    {
      return false;
    }
  else
    {
      return true;
    }
}

bool dir_exists(const std::string fname)
{
  return exists(fname);
}

bool file_exists(const std::string fname)
{
  return exists(fname);
}

bool regular_file(const std::string fname)
{
  fs::path p(fname);
  if(! fs::is_regular_file(p))
    {
      return false;
    }
  else
    {
      return true;
    }
}

bool file_can_read(const std::string fname)
{
  std::ifstream test(fname.c_str(), std::ios::in);
  bool ret;
  if(! test)
    {
      ret = false;
    }
  else
    {
      ret = true;
    }
  test.close();
  return ret;
}

bool file_can_create(const std::string fname)
{

  std::ofstream test(fname.c_str(), std::ios::out);
  bool ret;
  if(! test)
    {
      ret = false;
    }
  else
    {
      ret = true;
    }
  test.close();
  return ret;
}

bool dir_can_write(const std::string pname)
{
  fs::path p(pname);
  fs::path test = p;
  bool ret;
  test /= ".seidr_write_test";
  unsigned long counter = 0;
  while( fs::exists(test) )
    {
      test = p;
      test /= (".seidr_write_test" + std::to_string(counter));
      counter++;
    }
  std::ofstream wtest(test.string().c_str(), std::ios::out);
  if(! wtest)
    {
      ret = false;
    }
  else
    {
      ret = true;
      fs::remove(test);
    }
  wtest.close();
  return ret;
}

std::string to_absolute(const std::string xname)
{
  std::string ret = fs::absolute(xname).string();
  return ret;
}

std::string to_canonical(const std::string xname)
{
  std::string ret = fs::canonical(xname).string();
  return ret;
}

std::string dirname(const std::string xname)
{
  fs::path p(xname);
  fs::path q( p.parent_path() );
  return q.string();
}

bool create_directory(std::string base, std::string extend)
{
  fs::path p(base);
  p /= extend;
  return fs::create_directory(p);
}

bool starts_with(std::string fname, std::string pattern)
{
  fs::path p(fname);
  std::string s = p.filename().string();
  std::string ss = s.substr(0,pattern.size());
  return ss.compare(pattern) == 0;
}

void rename(std::string lhs, std::string rhs)
{
  fs::rename(lhs, rhs);
}

void remove(std::string fname, bool recursive)
{
  if (recursive)
    fs::remove_all(fname);
  else
    fs::remove(fname);
} 
