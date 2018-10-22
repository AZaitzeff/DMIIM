% script to make mex-file for redistC
mex COMPFLAGS='$COMPFLAGS -O2 -Wall -std=c++11 -fPIC' -O redistz.cpp redist.cpp idarray2d.cpp heap.cpp toolbox.cpp
