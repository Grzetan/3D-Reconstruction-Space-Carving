#include <vector>
#include <array>
#include <map>

#define SCENE_SIZE 100
#define VOXEL_SIZE 0.2

typedef std::vector<bool> Voxels; // false = visible, true = notvisible (to avoid filling vector with default values)
typedef std::array<double, 3> Point;
typedef std::vector<Point> Vertices;
typedef std::vector<std::vector<size_t>> Faces;

enum CubeFaces{
    LEFT=0, RIGHT=1, TOP=2, BOTTOM=3, FRONT=4, BACK=5
};

typedef std::map<CubeFaces, bool> FaceVisMap;
typedef std::map<CubeFaces, bool>::iterator FaceVisMapIterator;

const std::tuple<int, int, int> LEFT_FACE = {-1, 0, 0};
const std::tuple<int, int, int> RIGHT_FACE = {1, 0, 0};
const std::tuple<int, int, int> TOP_FACE = {0, 1, 0};
const std::tuple<int, int, int> BOTTOM_FACE = {0, -1, 0};
const std::tuple<int, int, int> FRONT_FACE = {0, 0, -1};
const std::tuple<int, int, int> BACK_FACE = {0, 0, 1};

const std::vector<std::tuple<int, int, int>> GUIDE = {LEFT_FACE, RIGHT_FACE, TOP_FACE, BOTTOM_FACE, FRONT_FACE, BACK_FACE};