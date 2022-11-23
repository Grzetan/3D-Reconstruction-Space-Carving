#pragma once

#include <iostream>
#include <map>
#include <vector>
#include "types.h"
#include "happly.h"
#include "bmplib.h"
#include <cmath>

void generateVoxels(Voxels& voxels);

void generateVoxelsMarchingCubes(Voxels& voxels);

size_t vIdx(Key key, UniqueVertices& uniqueVertices);

bool isCubeBackground(int x, int y, int z, Voxels& voxels);

int idx(int x, int y, int z, Voxels& voxels);

bool getVoxelValue(Voxels& voxels, int x, int y, int z);

bool rayAABBIntersection(const Ray& ray,
                         const AABB& box,
                         double& tstart,
                         double& tend);

void rayGridTraversal(Ray& ray, Voxels& voxels, const AABB& box);

void rotateAroundZAxis(Vec3& v, double angle);

Vec3 rotateUsingQuaterion(Vec3& v, double angle, RotationType type);

double degrees2radians(double degrees);

Bbox getBoundingBox(Voxels& voxels);

Cylinder getCylinder(Voxels& voxels, const Bbox& bbox, const Vec3& cameraPos);

void removeSingleVoxels(Voxels& voxels);

void createGroupVoxels(Voxels& voxels, int x, int y, int z, std::vector<VoxelArea>& groups);

void createGroupImage(BMP& image, int x, int y, std::vector<ImageArea>& groups);

void segmentImage(BMP& image);

void removeGate(Voxels& voxels, std::vector<std::array<int, 6>>& areasToRemove);

bool isPixelBackground(Pixel& p);

void calibrateVoxels(BMP& img, Voxels& voxels);