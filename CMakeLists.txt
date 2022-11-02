cmake_minimum_required(VERSION 3.0)

project(RabbitKSSD)

set(CMAKE_INSTALL_PREFIX ..)
option(USE_DEBUG "print the debug information" ON)
option(USE_Timer "print the time information" ON)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
#message("openmp found")

set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set(EXECUTABLE_OUTPUT_PATH .)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-g -O3 -march=native -lz ${CMAKE_CXX_FLAGS}")


if(USE_DEBUG)
	set(CMAKE_CXX_FLAGS "-D DEBUG ${CMAKE_CXX_FLAGS}")
endif()

if(USE_Timer)
	set(CMAKE_CXX_FLAGS "-D Timer ${CMAKE_CXX_FLAGS}")
endif()

aux_source_directory(src DIR_SRCS)
add_executable(kssd ${DIR_SRCS})

install(TARGETS kssd DESTINATION ./)






