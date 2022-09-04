#include <iostream>
#include "happly.h"
#include "types.h"
#include "params.h"
#include "utils.h"
#include <map>
#include <vector>

int main(int argc, char *argv[]){
    Params p(argc, argv);

    Voxels voxels = {};
    voxels.reserve(SCENE_SIZE*SCENE_SIZE*SCENE_SIZE);
    voxels[0] = true;

    // generateVoxels(voxels);

    Ray ray({0.5,0,0}, {-1,0,0.5});
    AABB box({0,0,0}, {VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE});

    double start, end;
    bool i = rayAABBIntersection(ray, box, start, end);

    std::cout << i << ", start: " << start << ", end: " << end << std::endl;
    return 0;
}