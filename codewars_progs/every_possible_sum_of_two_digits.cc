#include <string>
#include <vector>
#include "all_codewars.h"

std::vector<int> codewars_mysol::digits(int n) {
  const std::string word = std::to_string(n);
  std::vector<size_t> digits;
  for (const char& c : word)
  {
    digits.push_back((size_t)c - size_t('0'));
  }
  std::vector<int> arr_two_digits_sum;
  for (size_t i = 0; i != digits.size() - 1; ++i)
  {
    for (size_t j = i + 1; j != digits.size(); ++j)
    {
      arr_two_digits_sum.push_back((int)(digits[i] + digits[j]));
    }
  }
  return arr_two_digits_sum;
}
