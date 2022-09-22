#ifndef SHARED_LIB_H
#define SHARED_LIB_H

#ifdef __cplusplus
    extern "C"{
#endif

#ifdef BUILD_MROLLER_GENERATOR_DLL
    #define SHARED_LIB __declspec(dllexport)
#else
    #define SHARED_LIB __declspec(dllimport)
#endif

//declaration of exported methods

#ifdef __cplusplus
    }
#endif

#endif //end of DLL