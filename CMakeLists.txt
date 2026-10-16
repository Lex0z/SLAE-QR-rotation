cmake_minimum_required(VERSION 3.23)
project(lab2 C)

set(CMAKE_C_STANDARD 11)

add_executable(lab2 main.c
        Matrixs.c
        Matrixs.h
        LU-factorisation.c
        LU-factorisation.h)
