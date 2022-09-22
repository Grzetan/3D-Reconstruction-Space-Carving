g++ -c -DSPACE_CARVING_EXPORTS shared_lib.cpp -lstdc++fs -std=c++17
g++ -shared -o -shared_lib.dll shared_lib.o -Wl,--out-implib,libshared_lib.a -I/src