cmake_minimum_required(VERSION 3.5)
project(sualbsp)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS "-DIL_STD -g -ggdb -Wall -Wno-unknown-pragmas -fexceptions -std=c++20")

find_package(fmt REQUIRED)
find_package(cplex REQUIRED)
find_package(Boost CONFIG 1.80.0 REQUIRED COMPONENTS program_options timer system chrono filesystem)

include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
include_directories(AFTER .)
link_directories(${Boost_LIBRARY_DIRS})

set(EXTRA_LIBRARIES rt stdc++fs dl)

set(support_code random.cpp timer.cpp options.cpp util.cpp logging.cpp)

add_executable(hs hs.cpp instance.cpp hoffmann.cpp prbp.cpp solution.cpp dw.cpp lap.cpp sop.cpp models.cpp bbr.cpp lowerbounds.cpp ${support_code})
target_link_libraries(hs ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_LIBRARIES})
