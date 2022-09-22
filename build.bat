g++ -c -DBUILD_MROLLER_GENERATOR_DLL shared_lib.cpp
g++ -shared -o -shared_lib.dll shared_lib.o -Wl, --out-implib, libshared_lib.a