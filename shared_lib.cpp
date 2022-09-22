#include "shared_lib.h"

static double _radius;
static double _height;

void createBoundingCylinder(const char* path, unsigned int sceneSize, double rlVoxelSize, double xFOV, double yFOV)
{
    OutCylinder cylinder = spaceCarve(path, sceneSize, rlVoxelSize, xFOV, yFOV);
}

double getCylinderRadius()
{
    return _radius;
}

double getCylinderHeight()
{
    return _height;
}