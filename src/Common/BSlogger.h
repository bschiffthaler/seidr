/* * Seidr - Create and operate on gene crowd networks
 * Copyright (C) 2016-2019 Bastian Schiffthaler <b.schiffthaler@gmail.com>
 *
 * This file is part of Seidr.
 *
 * Seidr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Seidr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Seidr.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <iostream>
#include <ctime>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#define LOG_SILENT 0
#define LOG_ERR 1
#define LOG_WARN 2
#define LOG_INFO 3
#define LOG_TIME 4
#define LOG_DEBUG 5

#ifndef DEBUG
#define LOG_DEFAULT 4
#else
#define LOG_DEFAULT 5
#endif

#define LOG_INIT_COUT() logger log(std::cout, __PRETTY_FUNCTION__)
#define LOG_INIT_CERR() logger log(std::cerr, __PRETTY_FUNCTION__)
#define LOG_INIT_CLOG() logger log(std::clog, __PRETTY_FUNCTION__)
#define LOG_INIT_CUSTOM(X) logger log( (X), __PRETTY_FUNCTION__)

#ifdef BSLOG_NO_COLORS

#define BSLOG_TIME    "[ TIME    ]"
#define BSLOG_DEBUG   "[ DEBUG   ]"
#define BSLOG_ERROR   "[ ERROR   ]"
#define BSLOG_WARNING "[ WARNING ]"
#define BSLOG_INFO    "[ INFO    ]"

#else

#define BSLOG_TIME    "\033[0;35m[ TIME    ]\033[0;0m"
#define BSLOG_DEBUG   "[ DEBUG   ]"
#define BSLOG_ERROR   "\033[0;31m[ ERROR   ]\033[0;0m"
#define BSLOG_WARNING "\033[0;33m[ WARNING ]\033[0;0m"
#define BSLOG_INFO    "\033[0;34m[ INFO    ]\033[0;0m"

#endif

class logger  {
 public:
  logger(std::ostream&, unsigned, std::string);
  logger(std::ostream&, std::string n);
  template<typename T>
    friend logger& operator<<(logger& l, const T& s);
  logger& operator()(unsigned ll);
  void add_snapshot(std::string n, bool quiet = true) {
    time_t now; time(&now); _snaps.push_back(now);
    _snap_ns.push_back(n);
    if(_loglevel >= LOG_TIME && ! quiet)
      _fac << BSLOG_TIME << prep_time(*this) <<
	prep_name(*this) << ": Added snap '" << n << "'\n";
  }
  void set_log_level(unsigned ll) { _loglevel = ll;}
  void flush() { _fac.flush(); }
  void time_since_start();
  void time_since_last_snap();
  void time_since_snap(const std::string&);
  friend std::string prep_level(logger& l);
  friend std::string prep_time(logger& l);
  friend std::string prep_name(logger& l);
  static unsigned _loglevel;
 private:
  time_t _now;
  time_t _start;
  std::vector< time_t > _snaps;
  std::vector< std::string > _snap_ns;
  unsigned _message_level;
  std::ostream& _fac;
  std::string _name;
};

std::string prep_level(logger& l);
std::string prep_time(logger& l);
std::string prep_name(logger& l);


template<typename T>
logger& operator<<(logger& l, const T& s)
{
  if(l._message_level <= l._loglevel )
    {
      l._fac << s;
      return l;
    }
  else
    {
      return l;
    }
}
