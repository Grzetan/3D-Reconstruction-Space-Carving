#ifndef SHARED_LIB_H
#define SHARED_LIB_H

#ifdef __cplusplus
    extern "C"{
#endif

#include "src/api.h"

#ifdef BUILD_SPACE_CARVING_DLL
    #define SHARED_LIB __declspec(dllexport)
#else
    #define SHARED_LIB __declspec(dllimport)
#endif

//declaration of exported methods

/// @brief Space carving
/// @param path Path to folder with images (images should be .bmp and named incremental ex 0.bmp, 1.bmp, 2.bmp ...)
/// @param sceneSize <10, inf) Number of voxels in each dimention (The greater it is the slower and more accurate the algorithm is). Default value = 100
/// @param rlVoxelSize (0, inf) Size of one voxel in real world
/// @param xFOV <1, 180> Camera's field od view in X axis
/// @param yFOV <1, 180> Camera's field of view in Y axis
/// @return Cylinder class with two properties: height and radius
OutCylinder createBoundingCylinder(const char* path, unsigned int sceneSize, double rlVoxelSize, double xFOV, double yFOV);

/// @brief Space carving
/// @param path Path to folder with images (images should be .bmp and named incremental ex 0.bmp, 1.bmp, 2.bmp ...)
/// @return Cylinder class with two properties: height and radius
OutCylinder createBoundingCylinder(const char* path);

#ifdef __cplusplus
    }
#endif

#endif //end of DLL