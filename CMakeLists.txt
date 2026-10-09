cmake_minimum_required(VERSION 3.16)
project(EllipticFitting LANGUAGES CXX)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(MPI REQUIRED COMPONENTS CXX)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Boost::boost Eigen3::Eigen MPI::MPI_CXX)

add_executable(main_acl main_acl.cpp)
target_link_libraries(main_acl PRIVATE Boost::boost)

foreach(target main main_acl)
  set_target_properties(${target} PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
  )
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -pedantic)
  endif()
endforeach()
