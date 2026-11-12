#pragma once
#include <cmath>
#include <stdexcept>
#include <string>

inline long double parse_real(const char *text)
{
  std::size_t end = 0;
  const std::string value(text);
  const long double result = std::stold(value, &end);
  if (end != value.size() || !std::isfinite(result))
    throw std::invalid_argument("Expected a finite number: " + value);
  return result;
}

inline int parse_LCN(const char *text)
{
  const std::string value(text);
  if (value != "4" && value != "6")
    throw std::invalid_argument("LCN must be 4 or 6");
  return value == "4" ? 4 : 6;
}

inline void validate_k(long double k)
{
  if (!(0 < k && k < 1))
    throw std::invalid_argument("Expected 0 < k < 1; try different initial values");
}
