#include "shared_lib.h"

OutCylinder createBoundingCylinder(const char* path)
{
    return spaceCarve(path);
}


OutCylinder createBoundingCylinder(const char* path, unsigned int sceneSize, double rlVoxelSize, double xFOV, double yFOV)
{
    return spaceCarve(path, sceneSize, rlVoxelSize, xFOV, yFOV);
}