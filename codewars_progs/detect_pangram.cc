#include <string>
#include <iostream>
#include "all_codewars.h"

bool codewars_mysol::is_pangram(const std::string& s) {
    std::vector<bool> arr(26);
    for (const char& c : s)
    {
      const int c1 = (int)c;
      if (c1 >= 65 && c1 <= 90)
      {
        arr[c1 - 65] = true;
      } else if (c1 >= 97 && c1 <= 122)
      {
        arr[c1 - 97] = true;
      }
    }
    bool check = true;
    for (const bool& val : arr)
    {
      check = check && val;
    }
    return check;
}
