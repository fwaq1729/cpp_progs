#include <string>
#include <cctype>
#include "all_codewars.h"

std::string codewars_mysol::to_camel_case(std::string text) {
  std::string out = "";
  std::string ac = "";
  size_t count = 0;
  for (const char& c : text)
  {
     if ((c == '_') || (c == '-'))
     {
       if (count > 0)
       {
         ac[0] = std::toupper(ac[0]);
       }
       out += ac;
       ac = "";
       ++count;
     } else
     {
        ac += c;
     }
  }
  ac[0] = std::toupper(ac[0]);
  out += ac;
  return out;
}
