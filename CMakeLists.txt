cmake_minimum_required(VERSION 3.28)
project(MyRootProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Usa root-config per ottenere flag e librerie
execute_process(
    COMMAND root-config --cflags
    OUTPUT_VARIABLE ROOT_CXX_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND root-config --libs
    OUTPUT_VARIABLE ROOT_LIBRARIES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND root-config --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Aggiungi include e flag
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Crea l'eseguibile
add_executable(prima_prova main.cpp)

# imposta opzioni di compilazione (separate gli argomenti)
string(REPLACE " " ";" ROOT_CFLAGS_LIST "${ROOT_CXX_FLAGS}")
target_compile_options(prima_prova PRIVATE ${ROOT_CFLAGS_LIST})

# Linka con le librerie ROOT
target_link_libraries(prima_prova ${ROOT_LIBRARIES})
# --- Aggiungi un target "run" per eseguire il programma ---
add_custom_target(run
    COMMAND prima_prova
    DEPENDS prima_prova
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Esecuzione del programma ROOT..."
)

