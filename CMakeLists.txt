cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(signal)

if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

add_executable(signal signal.C)
target_link_libraries(signal Garfield::Garfield)
