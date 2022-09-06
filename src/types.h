#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <cmath>

#define SCENE_SIZE 100
#define VOXEL_SIZE 0.2

class Vec3{
public:
    double x, y, z;

    Vec3(){
        x = 0;
        y = 0;
        z = 0;
    }

    Vec3(double x_, double y_, double z_){
        x = x_;
        y = y_;
        z = z_;
    }

    Vec3 operator*(double t){
        return {x*t, y*t, z*t};
    }

    Vec3 operator+(Vec3 v){
        return {x+v.x, y+v.y, z+v.z};
    }

    Vec3 operator-(Vec3 v){
        return {x-v.x, y-v.y, z-v.z};
    }

    Vec3 operator+(double a){
        return {x+a, y+a, z+a};
    }

    double len(){
        return std::sqrt(x*x + y*y + z*z);
    }

    void normalize(){
        double length = len();
        x /= length;
        y /= length;
        z /= length;
    }
};

typedef std::vector<bool> Voxels; // false = visible, true = notvisible (to avoid filling vector with default values)
typedef std::vector<std::array<double, 3>> Vertices;
typedef std::vector<std::vector<size_t>> Faces;

enum CubeFaces{
    LEFT=0, RIGHT=1, TOP=2, BOTTOM=3, FRONT=4, BACK=5
};

typedef std::map<CubeFaces, bool> FaceVisMap;
typedef std::map<CubeFaces, bool>::iterator FaceVisMapIterator;

typedef std::tuple<double, double, double> Key;
typedef std::map<Key, int> UniqueVertices;

const std::tuple<int, int, int> LEFT_FACE = {-1, 0, 0};
const std::tuple<int, int, int> RIGHT_FACE = {1, 0, 0};
const std::tuple<int, int, int> TOP_FACE = {0, 1, 0};
const std::tuple<int, int, int> BOTTOM_FACE = {0, -1, 0};
const std::tuple<int, int, int> FRONT_FACE = {0, 0, -1};
const std::tuple<int, int, int> BACK_FACE = {0, 0, 1};

const std::vector<std::tuple<int, int, int>> GUIDE = {LEFT_FACE, RIGHT_FACE, TOP_FACE, BOTTOM_FACE, FRONT_FACE, BACK_FACE};

class Ray{
public:
    Vec3 dir, orig, invDir;
    int sign[3];

    Ray(Vec3 dir_, Vec3 orig_){
        dir = dir_;
        orig = orig_;
        invDir.x = 1.0 / dir.x;
        invDir.y = 1.0 / dir.y;
        invDir.z = 1.0 / dir.z;
        sign[0] = invDir.x < 0;
        sign[1] = invDir.y < 0;
        sign[2] = invDir.z < 0;
    }
};

class AABB{
public:
    Vec3 bounds[2];
    AABB(Vec3 min, Vec3 max){
        bounds[0] = min;
        bounds[1] = max;
    }
};