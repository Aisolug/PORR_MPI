cmake_minimum_required(VERSION 3.9)
project(porr2)

set(CMAKE_CXX_STANDARD 11)


add_executable(porr2
        jacobi.c
        jacobi.h
        main.c
        richardson.c
        richardson.h
        testFun.c
        testFun.h)

target_link_libraries(porr2 m)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(porr2 ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(porr2 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(porr2 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()