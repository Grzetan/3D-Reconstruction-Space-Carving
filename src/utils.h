#pragma once

#include <iostream>
#include <map>
#include <vector>
#include "types.h"
#include "happly.h"
#include <cmath>

void generateVoxels(Voxels& voxels);

size_t vIdx(Key key, UniqueVertices& uniqueVertices);

bool isCubeBackground(int x, int y, int z, Voxels& voxels);

int idx(int x, int y, int z);

bool rayAABBIntersection(const Ray& ray,
                         const AABB& box,
                         double& tstart,
                         double& tend);

void rayGridTraversal(Ray& ray, Voxels& voxels, const AABB& box);