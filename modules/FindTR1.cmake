# Check availability of unordered maps
cmake_policy(PUSH)
cmake_minimum_required(VERSION 2.6.3)
cmake_policy(POP)

INCLUDE (CheckCXXSourceCompiles)

check_cxx_source_compiles(
  "
  #include <tr1/unordered_map>
  int main() 
   {
     std::tr1::unordered_map<int, int> m;
     return 0;
   }
  "
TR1_UNORDERED_MAP_TR1_TR1)


check_cxx_source_compiles(
  "
  #include <tr1/unordered_map>
  int main() 
   {
     std::unordered_map<int, int> m;
     return 0;
   }
  "
TR1_UNORDERED_MAP_TR1_STD)

check_cxx_source_compiles(
  "
  #include <unordered_map>
  int main() 
   {
     std::tr1::unordered_map<int, int> m;
     return 0;
   }
  "
TR1_UNORDERED_MAP_STD_TR1)

check_cxx_source_compiles(
  "
  #include <unordered_map>
  int main() 
   {
     std::unordered_map<int, int> m;
     return 0;
   }
  "
TR1_UNORDERED_MAP_STD_STD)

if (TR1_UNORDERED_MAP_STD_STD)
  set (TR1_UNORDERED_MAP_FOUND true)
  set (TR1_HEADER_PREFIX false)
  set (TR1_NAMESPACE_PREFIX false)
elseif(TR1_UNORDERED_MAP_TR1_TR1)
  set (TR1_UNORDERED_MAP_FOUND true)
  set (TR1_HEADER_PREFIX true)
  set (TR1_NAMESPACE_PREFIX true)
elseif(TR1_UNORDERED_MAP_TR1_STD)
  set (TR1_UNORDERED_MAP_FOUND true)
  set (TR1_HEADER_PREFIX true)
  set (TR1_NAMESPACE_PREFIX false)
elseif(TR1_UNORDERED_MAP_STD_TR1)
  set (TR1_UNORDERED_MAP_FOUND true)
  set (TR1_HEADER_PREFIX false)
  set (TR1_NAMESPACE_PREFIX true)
endif()

