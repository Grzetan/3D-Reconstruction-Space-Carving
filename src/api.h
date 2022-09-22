#pragma once
#include "bmplib.h"
#include "types.h"
#include "utils.h"
#include <filesystem>
#include <algorithm>

OutCylinder spaceCarve(const char* path, unsigned int sceneSize, double rlVoxelSize);