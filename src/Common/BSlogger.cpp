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

#include <BSlogger.h>

#include <algorithm>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

unsigned logger::_loglevel = LOG_DEFAULT;

logger::logger(std::ostream& f, std::string n) :
  _now(0),
  _start(0),
  _message_level(LOG_SILENT),
  _fac(f),
  _name(std::move(n))
{
  time(&_now);
  time(&_start);
}

logger::logger(std::ostream& f, unsigned ll, std::string n) :
  _now(0),
  _start(0),
  _message_level(LOG_SILENT),
  _fac(f),
  _name(std::move(n))
{
  _loglevel = ll;
  time(&_now);
  time(&_start);
}

logger& logger::operator()(unsigned ll) {
  _message_level = ll;
  if (_message_level <= _loglevel )
  {
    _fac << prep_level(*this) << prep_time(*this) <<
         prep_name(*this) << ": ";
  }
  return *this;
}

std::string prep_level(logger& l)
{
  switch (l._message_level)
  {
  case LOG_ERR:
    return BSLOG_ERROR; break;
  case LOG_WARN:
    return BSLOG_WARNING; break;
  case LOG_INFO:
    return BSLOG_INFO; break;
  case LOG_DEBUG:
    return BSLOG_DEBUG; break;
  case LOG_TIME:
    return BSLOG_TIME; break;
  default:
    return "";
  }
  return "";
}

std::string prep_time(logger& l)
{
  time(&l._now);
  struct tm * t;
  t = localtime(&l._now);
  std::string s, m, h, D, M, Y;
  s = std::to_string(t->tm_sec);
  m = std::to_string(t->tm_min);
  h = std::to_string(t->tm_hour);
  D = std::to_string(t->tm_mday);
  M = std::to_string(t->tm_mon + 1);
  Y = std::to_string(t->tm_year + 1900);

  if (t->tm_sec < 10)
  {
    s = "0" + s;
  }
  if (t->tm_min < 10)
  {
    m = "0" + m;
  }
  if (t->tm_hour < 10)
  {
    h = "0" + h;
  }
  if (t->tm_mday < 10)
  {
    D = "0" + D;
  }
  if (t->tm_mon + 1 < 10)
  {
    M = "0" + M;
  }

  std::string ret = "[ " + Y + "-" + M + "-" + D +
                    "T" + h + ":" + m + ":" + s + " ]";

  return ret;
}

std::string prep_name(logger& l)
{
  return "[ " + l._name + " ]";
}

void logger::time_since_start()
{
  if (_loglevel >= LOG_TIME)
  {
    time(&_now);
    _message_level = LOG_TIME;
    _fac << prep_level(*this) << prep_time(*this) <<
         prep_name(*this) << ": " <<
         difftime(_now, _start) << "s since instantiation\n";
  }
}

void logger::time_since_last_snap()
{
  if (_loglevel >= LOG_TIME)
  {
    time(&_now);
    _message_level = LOG_TIME;
    _fac << prep_level(*this) << prep_time(*this) <<
         prep_name(*this) << ": " <<
         difftime(_now, _snaps.back()) << "s since snap '" <<
         _snap_ns.back() << "'\n";
  }
}

void logger::time_since_snap(const std::string& s)
{
  if (_loglevel >= LOG_TIME)
  {
    time(&_now);
    auto it = find(_snap_ns.begin(), _snap_ns.end(), s);
    uint64_t dist = std::distance(_snap_ns.begin(), it);

    _message_level = LOG_TIME;
    _fac << prep_level(*this) << prep_time(*this) <<
         prep_name(*this) << ": " <<
         difftime(_now, _snaps[dist]) << "s since snap '" <<
         _snap_ns[dist] << "'\n";
  }
}
