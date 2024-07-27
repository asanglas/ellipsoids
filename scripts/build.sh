mkdir -p build
FLAGS=$(cat compile_flags.txt)

clang src/main.c $FLAGS -o build/ellipsoids
