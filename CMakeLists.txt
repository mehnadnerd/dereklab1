cmake_minimum_required(VERSION 3.21)
project(lab1a)

set(CMAKE_CXX_STANDARD 23)

add_executable(lab1a gen_matrix.c gen_matrix.h bmm.c my_malloc.c my_malloc.h)
