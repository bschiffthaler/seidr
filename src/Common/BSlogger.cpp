#include <BSlogger.h>
#include <iostream>
#include <ctime>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

unsigned logger::_loglevel = LOG_DEFAULT;

logger::logger(std::ostream& f, std::string n) :
	      _message_level(LOG_SILENT),
	      _fac(f),
	      _name(n)
{
  time(&_now);
  time(&_start);
}

logger::logger(std::ostream& f, unsigned ll, std::string n) :
	      _message_level(LOG_SILENT),
	      _fac(f),
	      _name(n)
{
  time(&_now);
  time(&_start);
  _loglevel = ll;
}

logger& logger::operator()(unsigned ll){
  _message_level = ll;
  if(_message_level <= _loglevel )
    {
      _fac << prep_level(*this) << prep_time(*this) <<
	prep_name(*this) << ": ";
    }
  return *this;
}

std::string prep_level(logger& l)
{
  switch(l._message_level)
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
  std::string s,m,h,D,M,Y;
  s = std::to_string(t->tm_sec);
  m = std::to_string(t->tm_min);
  h = std::to_string(t->tm_hour);
  D = std::to_string(t->tm_mday);
  M = std::to_string(t->tm_mon + 1);
  Y = std::to_string(t->tm_year + 1900);

  if(t->tm_sec < 10) s = "0" + s;
  if(t->tm_min < 10) m = "0" + m;
  if(t->tm_hour < 10) h = "0" + h;
  if(t->tm_mday < 10) D = "0" + D;
  if(t->tm_mon + 1 < 10) M = "0" + M;
  
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
  if(_loglevel >= LOG_TIME)
    {
      time(&_now);
      _message_level = LOG_TIME;
      _fac << prep_level(*this) << prep_time(*this) <<
	prep_name(*this) << ": " << 
	difftime(_now,_start) << "s since instantiation\n";
    }
}

void logger::time_since_last_snap()
{
  if(_loglevel >= LOG_TIME)
    {
      time(&_now);
      _message_level = LOG_TIME;
      _fac << prep_level(*this) << prep_time(*this) <<
	prep_name(*this) << ": " <<
	difftime(_now,_snaps.back()) << "s since snap '" <<
	_snap_ns.back() << "'\n";
    }
}

void logger::time_since_snap(std::string s)
{
  if(_loglevel >= LOG_TIME)
    {
      time(&_now);
      auto it = find(_snap_ns.begin(), _snap_ns.end(), s);
      unsigned long dist = std::distance(_snap_ns.begin(), it);
      
      _message_level = LOG_TIME;
      _fac << prep_level(*this) << prep_time(*this) <<
	prep_name(*this) << ": " <<
	difftime(_now,_snaps[dist]) << "s since snap '" <<
	_snap_ns[dist] << "'\n";
    }
}
