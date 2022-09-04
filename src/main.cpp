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

    Ray ray({1,0,0}, {-10,1.8,0.3});
    AABB box({0,0,0}, {VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE});

    rayGridTraversal(ray, voxels, box);
    return 0;
}