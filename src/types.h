#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <cmath>

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

struct Voxels{
    int SCENE_SIZE;
    double VOXEL_SIZE, VOXEL_RL_SIZE;
    std::vector<bool> data; // false = visible, true = notvisible (to avoid filling vector with default values)

    Voxels(int s_size, double v_size, double vrl_size){
        SCENE_SIZE = s_size;
        VOXEL_SIZE = v_size;
        VOXEL_RL_SIZE = vrl_size;
        data.reserve(SCENE_SIZE*SCENE_SIZE*SCENE_SIZE);
    }
};

typedef std::vector<std::array<double, 3>> Vertices;
typedef std::vector<std::vector<size_t>> Faces;

enum RotationType{
    UP_DOWN, LEFT_RIGHT
};

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

class Quaternion{
public:
    double a, b, c ,d;

    Quaternion(double a_, double b_, double c_, double d_){
        a = a_;
        b = b_;
        c = c_;
        d = d_;
    }

    Quaternion(double angle, Vec3& v){
        create(angle, v);
    }

    Quaternion(Vec3 v){
        a = 0;
        b = v.x;
        c = v.y;
        d = v.z;
    }

    void create(double angle, Vec3& v){
        double tmp = std::sin(angle / 2.0f);
        a = std::cos(angle / 2.0f);
        b = tmp*v.x;
        c = tmp*v.y;
        d = tmp*v.z;
    }

    Quaternion inverse(){
        return {a, -b, -c, -d};
    }

    Quaternion hamilton(Quaternion& q){
        return {
            a*q.a - b*q.b - c*q.c - d*q.d,
            a*q.b + b*q.a + c*q.d - d*q.c,
            a*q.c - b*q.d + c*q.a + d*q.b,
            a*q.d + b*q.c - c*q.b + d*q.a
        };
    }
};