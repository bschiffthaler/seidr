#pragma once

#include <tclap/CmdLine.h>
#include "tclap/DocBookOutput.h"
#include "tclap/ZshCompletionOutput.h"
#include <list>
#include <string>
#include <sstream>

class tclap_seidr_format : public TCLAP::StdOutput
{
 public:
  virtual void usage(TCLAP::CmdLineInterface& c)
  {
    std::list<TCLAP::Arg*> args = c.getArgList();
    std::cout << "Usage:\n  "
              << c.getProgramName()
              << ' ';
    unsigned int counter = 2;
    int words = args.size();
    for (auto& a : args)
      {
        std::string sid = a->shortID();
        counter += sid.size() + 1;
        if(counter > 80)
          {
            std::cout << "\n  ";
            counter = 2;
          }
        words--;
        std::cout << sid << (words > -1 ? ' ' : '\n');
      }
    std::cout << "Where:\n";
    for (auto& a : args)
      {
        std::cout
          << "  -" << a->getFlag()
          << ", --" << a->getName()
          << "\n    ";
        counter = 4;
        words = 0;
        for(auto& c : a->getDescription())
            if(c == ' ')
              words++;
        std::stringstream ss(a->getDescription());
        std::string word;
        while(std::getline(ss, word, ' '))
          {
            counter += word.size() + 1;
            if(counter > 80)
              {
                std::cout << "\n    ";
                counter = 4;
              }
            words--;
            std::cout << word << (words > -1 ? ' ' : '\n');
          }
        std::cout << '\n';
      }
    std::cout << c.getMessage();
  }
};
