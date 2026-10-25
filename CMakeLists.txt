cmake_minimum_required(VERSION 3.23)
project(SLAE-LU-method C)

set(CMAKE_C_STANDARD 11)

add_executable(SLAE-LU-method main.c
        Matrixs.c
        Matrixs.h
        LU-factorisation.c
        LU-factorisation.h)

# Передаем путь к корню проекта в код как макрос PROJECT_ROOT_DIR
target_compile_definitions(${PROJECT_NAME} PRIVATE
        PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)
