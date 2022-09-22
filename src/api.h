#pragma once
#include "bmplib.h"
#include "types.h"
#include "utils.h"
#include <filesystem>
#include <algorithm>

/// @brief Space carving
/// @param path Path to folder with images (images should be .bmp and named incremental ex 0.bmp, 1.bmp, 2.bmp ...)
/// @param sceneSize <10, inf) Number of voxels in each dimention (The greater it is the slower and more accurate the algorithm is). Default value = 100
/// @param rlVoxelSize (0, inf) Size of one voxel in real world
/// @param xFOV <1, 180> Camera's field od view in X axis
/// @param yFOV <1, 180> Camera's field of view in Y axis
/// @return Cylinder class with two properties: height and radius
OutCylinder spaceCarve(const char* path, unsigned int sceneSize = 100, double rlVoxelSize = 1, double xFOV = 46.7, double yFOV = 46.7);