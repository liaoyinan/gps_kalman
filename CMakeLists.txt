cmake_minimum_required(VERSION 3.17)
project(_ C)

set(CMAKE_C_STANDARD 99)

include_directories(.)

add_executable(gps
        test_gps.c gps_kalman.c gps_kalman.h)
