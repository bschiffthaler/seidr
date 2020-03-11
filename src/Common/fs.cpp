//
// Seidr - Create and operate on gene crowd networks
// Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
//
// This file is part of Seidr.
//
// Seidr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Seidr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
//

#include <fstream>
#include <boost/filesystem.hpp>
#include <fs.h>

namespace fs = boost::filesystem;

bool exists(const std::string& fname)
{
  fs::path p(fname);
  return fs::exists(p);
}

bool dir_exists(const std::string& fname)
{
  return exists(fname);
}

bool file_exists(const std::string& fname)
{
  return exists(fname);
}

void assert_exists(const std::string& fname)
{
  if (! exists(fname))
  {
    throw std::runtime_error("No such file or directory: " + fname);
  }
}

void assert_no_overwrite(const std::string& fname)
{
  if (exists(fname))
  {
    throw std::runtime_error("File exists: " + fname);
  }
}

void assert_can_read(const std::string& fname)
{
  if (! file_can_read(fname))
  {
    throw std::runtime_error("Cannot open file for reading: " + fname);
  }
}

void assert_dir_is_writeable(const std::string& fname)
{
  if (! dir_can_write(fname))
  {
    throw std::runtime_error("Cannot create files in directory: " + fname);
  }
}

void assert_is_regular_file(const std::string& fname)
{
  if (! regular_file(fname))
  {
    throw std::runtime_error("Not a regular file: " + fname);
  }
}

bool regular_file(const std::string& fname)
{
  fs::path p(fname);
  return fs::is_regular_file(p);
}

bool file_can_read(const std::string& fname)
{
  std::ifstream test(fname.c_str(), std::ios::binary);
  bool ret = test.good();
  test.close();
  return ret;
}

bool file_can_create(const std::string& fname)
{
  std::ofstream test(fname.c_str(), std::ios::out);
  bool ret = test.good();
  test.close();
  return ret;
}

bool dir_can_write(const std::string& pname)
{
  fs::path p(pname);
  fs::path test = p;
  test /= fs::unique_path();
  std::ofstream wtest(test.string().c_str(), std::ios::out);
  bool ret = wtest.good();
  wtest.close();
  if (exists(test.string()))
  {
    remove(test.string());
  }
  return ret;
}

std::string to_absolute(const std::string& xname)
{
  std::string ret = fs::absolute(xname).string();
  return ret;
}

std::string to_canonical(const std::string& xname)
{
  std::string ret = fs::canonical(xname).string();
  return ret;
}

std::string dirname(const std::string& xname)
{
  fs::path p(xname);
  fs::path q( p.parent_path() );
  return q.string();
}

bool create_directory(const std::string& path)
{
  fs::path p(path);
  return fs::create_directory(p);
}

bool create_directory(const std::string& base, const std::string& extend)
{
  fs::path p(base);
  p /= extend;
  return fs::create_directory(p);
}

bool starts_with(const std::string& fname, const std::string& pattern)
{
  fs::path p(fname);
  std::string s = p.filename().string();
  std::string ss = s.substr(0, pattern.size());
  return ss == pattern;
}

void rename(const std::string& lhs, const std::string& rhs)
{
  try
  {
    fs::rename(lhs, rhs);
  }
  // Fall back to copying in case relinking doesn't work
  catch (const fs::filesystem_error& e)
  {
    fs::copy_file(lhs, rhs, fs::copy_option::overwrite_if_exists);
    fs::remove(lhs);
  }
}

void remove(const std::string& fname, bool recursive)
{
  if (recursive)
  {
    fs::remove_all(fname);
  }
  else
  {
    fs::remove(fname);
  }
}

std::string replace_ext(const std::string& fname, const std::string& new_ext)
{
  fs::path p(fname);
  p.replace_extension(new_ext);
  return p.string();
}

std::string tempfile(const std::string& tempdir)
{
  fs::path tempf;
  boost::filesystem::path dir;
  dir = fs::canonical(dir);

  if (tempdir.empty())
  {
    dir = fs::temp_directory_path();
  }
  else
  {
    dir = fs::path(tempdir);
  }
  
  // Guarantee that we are not overwriting anything
  do
  {
    auto tfile = fs::unique_path();
    tempf = (dir /= tfile);
  }
  while(fs::exists(tempf));


  return tempf.string();
}

std::string tempdir(const std::string& tempdir)
{
  boost::filesystem::path dir;
  if (tempdir.empty())
  {
    dir = fs::temp_directory_path();
  }
  else
  {
    dir = fs::path(tempdir);
  }
  dir = fs::canonical(dir);
  return dir.string();
}
