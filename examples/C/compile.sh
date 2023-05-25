g++ --std=c++17 -O3 -I../../include/ -I../../external/fmt/include/ -c ../../src/SHARPlib/*.cpp ../../src/SHARPlib/CWrap/*.cpp
ar rcs SHARPlib_C.a *.o
gcc -O3 -I../../include -o main main.c SHARPlib_C.a -lstdc++ -lm
