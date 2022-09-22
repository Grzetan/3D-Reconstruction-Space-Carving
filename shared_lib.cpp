#include "shared_lib.h"

OutCylinder createBoundingCylinder(const char* path)
{
    unsigned int sceneSize = 100;
    double rlVoxelSize = 1;
    double xFOV = 46.7;
     double yFOV = 46.7;
    return spaceCarve(path, sceneSize, rlVoxelSize, xFOV, yFOV);
}


OutCylinder createBoundingCylinder(const char* path, unsigned int sceneSize, double rlVoxelSize, double xFOV, double yFOV)
{
    return spaceCarve(path, sceneSize, rlVoxelSize, xFOV, yFOV);
}