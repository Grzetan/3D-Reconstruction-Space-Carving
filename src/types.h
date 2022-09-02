#include <vector>
#include <array>

#define SCENE_SIZE 100
#define VOXEL_SIZE 0.2

typedef std::vector<bool> Voxels;
typedef std::array<double, 3> Point;
typedef std::vector<Point> Vertices;
typedef std::vector<std::vector<size_t>> Faces;