cmake_minimum_required(VERSION 3.5)

# set the project name
project(BFS VERSION 0.1)

include_directories("${PROJECT_BINARY_DIR}")

if(MSVC)
  add_compile_options(/W4 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_library(mcl_DIR NAMES mcl REQUIRED)
add_library(mcl SHARED IMPORTED)
set_target_properties(mcl PROPERTIES
  IMPORTED_LOCATION "${mcl_DIR}"
)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
link_libraries(Threads::Threads)

add_library(utils OBJECT utils.cpp)
add_library(graph OBJECT Graph.cpp)
add_library(result_database OBJECT ResultDatabase.cpp)
add_library(option_parser OBJECT OptionParser.cpp Option.cpp)

add_custom_target(
        bfs_kernel ALL
        ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/bfs_iiit.cl
            ${CMAKE_CURRENT_BINARY_DIR}/bfs_iiit.cl )
                
# add the executables
add_executable(bfs main.cpp BFS.cpp
                $<TARGET_OBJECTS:utils> $<TARGET_OBJECTS:graph> 
                $<TARGET_OBJECTS:result_database> $<TARGET_OBJECTS:option_parser>)

# add_executable(bfs_modified  main.cpp BFS_modified.cpp
#                 $<TARGET_OBJECTS:utils> $<TARGET_OBJECTS:graph> 
#                 $<TARGET_OBJECTS:result_database> $<TARGET_OBJECTS:option_parser>)

target_link_libraries(bfs mcl OpenCL m rt)
#target_link_libraries(bfs_modified mcl OpenCL m rt)