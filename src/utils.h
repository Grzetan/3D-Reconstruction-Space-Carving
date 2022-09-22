#pragma once

#include <iostream>
#include <map>
#include <vector>
#include "types.h"
#include "happly.h"
#include <cmath>

void generateVoxels(Voxels& voxels);

void generateVoxelsMarchingCubes(Voxels& voxels);

size_t vIdx(Key key, UniqueVertices& uniqueVertices);

bool isCubeBackground(int x, int y, int z, Voxels& voxels);

int idx(int x, int y, int z, Voxels& voxels);

bool rayAABBIntersection(const Ray& ray,
                         const AABB& box,
                         double& tstart,
                         double& tend);

void rayGridTraversal(Ray& ray, Voxels& voxels, const AABB& box);

void rotateAroundZAxis(Vec3& v, double angle);

Vec3 rotateUsingQuaterion(Vec3& v, double angle, RotationType type);

double degrees2radians(double degrees);

Vec3 getBoundingBox(Voxels& voxels);

Cylinder getCylinder(Voxels& voxels, Vec3& bbox);

OutCylinder getOutCylinder(Voxels& voxels, Vec3& bbox);