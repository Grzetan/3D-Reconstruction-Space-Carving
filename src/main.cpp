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

    Ray ray({1,0,0}, {-10,0.5,0.3});
    AABB box({0,0,0}, {VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE});

    rayGridTraversal(ray, voxels, box);
    std::cout << "DONE Traversing" << std::endl;

    generateVoxels(voxels);
    return 0;
}