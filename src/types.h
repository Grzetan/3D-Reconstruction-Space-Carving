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

    double dist(Vec3& b){
        return std::sqrt(std::pow(b.x - x, 2) + std::pow(b.y - y, 2) + std::pow(b.z - z, 2));
    }

    void print(){
        std::cout << x << ", " << y << ", " << z << std::endl;
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

struct Vec3Int{
    int x,y,z;
};

struct VoxelArea{
    Vec3Int start;
    Vec3Int end;

    bool isValid(int minGroupSize = 7){
        return end.x - start.x > minGroupSize || end.y - start.y > minGroupSize || end.z - start.z > minGroupSize;
    }
};

struct ImageArea{
    std::array<int, 2> start;
    std::array<int, 2> end;

    bool isValid(int minGroupSize = 50){
        return end[0] - start[0] > minGroupSize || end[1] - start[1] > minGroupSize;
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
typedef std::tuple<int, int, int> FaceKey; 

const FaceKey LEFT_FACE = {-1, 0, 0};
const FaceKey RIGHT_FACE = {1, 0, 0};
const FaceKey TOP_FACE = {0, 1, 0};
const FaceKey BOTTOM_FACE = {0, -1, 0};
const FaceKey FRONT_FACE = {0, 0, -1};
const FaceKey BACK_FACE = {0, 0, 1};

const std::vector<FaceKey> GUIDE = {LEFT_FACE, RIGHT_FACE, TOP_FACE, BOTTOM_FACE, FRONT_FACE, BACK_FACE};

const std::vector<std::array<int, 4>> ACTIVE_VERTICES = {
    {0,3,7,4}, // Left
    {6,5,2,1}, // Right
    {7,6,5,4}, // Top
    {3,0,1,2}, // Bottom
    {7,6,2,3}, // Front
    {4,0,1,5}  // Back
};

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

struct Cylinder{
    double r;
    double h;
    Vec3 center;
};

struct OutCylinder{
    double radius;
    double height;
};

struct Bbox{
    Vec3 min;
    Vec3 max;
};

// Relative coordinates of middle point on given edge
const float edgeTable[12][3] = {
    {0.5, 0.0, 1.0}, // 0
    {1.0, 0.0, 0.5}, // 1
    {0.5, 0.0, 0.0}, // 2
    {0.0, 0.0, 0.5}, // 3
    {0.5, 1.0, 1.0}, // 4
    {1.0, 1.0, 0.5}, // 5
    {0.5, 1.0, 0.0}, // 6
    {0.0, 1.0, 0.5}, // 7
    {0.0, 0.5, 1.0}, // 8
    {1.0, 0.5, 1.0}, // 9
    {1.0, 0.5, 0.0}, // 10
    {0.0, 0.5, 0.0}, // 11
};