cmake_minimum_required(VERSION 3.13)

# ! INSERT PROJECT NAME BELOW !

Set(This PROJECT_NAME)

Project(${This})

set(CMAKE_CXX_STANDARD_17)

enable_testing()

add_subdirectory(src)
add_subdirectory(lib/googletest)
add_subdirectory(tst)

include_directories(src)
