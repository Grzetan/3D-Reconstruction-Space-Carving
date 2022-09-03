#include <iostream>
#include "happly.h"
#include "types.h"
#include "params.h"
#include <map>
#include <vector>

int idx(int x, int y, int z){
    x = std::min(SCENE_SIZE-1, std::max(0, x));
    y = std::min(SCENE_SIZE-1, std::max(0, y));
    z = std::min(SCENE_SIZE-1, std::max(0, z));

    return x + y*SCENE_SIZE + z*SCENE_SIZE*SCENE_SIZE;
}

bool isCubeBackground(int x, int y, int z, Voxels& voxels){
    if(x<0 || x>SCENE_SIZE-1 || y<0 || y>SCENE_SIZE-1 || z<0 || z>SCENE_SIZE-1){
        return true;
    }

    return voxels[idx(x, y, z)];
};

// for each face check if it should be displayed:
// 

void generateVoxels(Voxels& voxels){
    happly::PLYData out;

    Vertices vertices;
    Faces faces;

    int x,y,z;
    int vis = 0;

    for(x=0; x<SCENE_SIZE; x++){
        for(y=0; y<SCENE_SIZE; y++){
            for(z=0; z<SCENE_SIZE; z++){
                // Draw voxel if val is true and there is at least one neighbour not visible
                if(voxels[idx(x, y, z)] == false){
                    // If voxel is on the edge, skip neightbour test
                    // Check if all of the neighbours are visible, if so skip
                    
                    FaceVisMap faceVisibility;
                    faceVisibility[LEFT] = false;
                    faceVisibility[RIGHT] = false;
                    faceVisibility[TOP] = false;
                    faceVisibility[BOTTOM] = false;
                    faceVisibility[FRONT] = false;
                    faceVisibility[BACK] = false;

                    for(int i=0; i<GUIDE.size(); i++){
                        int cx, cy, cz;
                        std::tie(cx, cy, cz) = GUIDE[i];
                        faceVisibility[CubeFaces(i)] = !isCubeBackground(x+cx, y+cy, z+cz, voxels);
                    }

                    FaceVisMapIterator it;
                    for(it = faceVisibility.begin(); it != faceVisibility.end(); it++){
                        if(!it->second) vis++;
                    }

                    // Front, bottom, left
                    double mainX = x * VOXEL_SIZE, mainY = y * VOXEL_SIZE, mainZ = z * VOXEL_SIZE;
                    double offsetX = mainX + VOXEL_SIZE, offsetY = mainY + VOXEL_SIZE, offsetZ = mainZ + VOXEL_SIZE;

                    const size_t relVerticesSize = vertices.size();

                    // Generate vertex only if neighbour face is visible

                    if(!faceVisibility[FRONT] || !faceVisibility[BOTTOM] || !faceVisibility[LEFT]){
                        vertices.push_back({mainX, mainY, mainZ}); // 0: Front bottom left
                    }

                    if(!faceVisibility[FRONT] || !faceVisibility[TOP] || !faceVisibility[RIGHT]){
                        vertices.push_back({mainX, offsetY, mainZ}); // 1: Front top left
                    }

                    if(!faceVisibility[FRONT] || !faceVisibility[BOTTOM] || !faceVisibility[RIGHT]){
                        vertices.push_back({offsetX, mainY, mainZ}); // 2: Front bottom right
                    }

                    if(!faceVisibility[FRONT] || !faceVisibility[TOP] || !faceVisibility[RIGHT]){
                        vertices.push_back({offsetX, offsetY, mainZ}); // 3: Front top right
                    }

                    if(!faceVisibility[BACK] || !faceVisibility[BOTTOM] || !faceVisibility[LEFT]){
                        vertices.push_back({mainX, mainY, offsetZ}); // 4: Back bottom left
                    }

                    if(!faceVisibility[BACK] || !faceVisibility[TOP] || !faceVisibility[LEFT]){
                        vertices.push_back({mainX, offsetY, offsetZ}); // 5: Back top left
                    }

                    if(!faceVisibility[BACK] || !faceVisibility[BOTTOM] || !faceVisibility[RIGHT]){
                        vertices.push_back({offsetX, mainY, offsetZ}); // 6: Back bottom right
                    }

                    if(!faceVisibility[BACK] || !faceVisibility[TOP] || !faceVisibility[RIGHT]){
                        vertices.push_back({offsetX, offsetY, offsetZ}); // 7: Back top right
                    }

                    // // Generate faces
                    // faces.push_back({relVerticesSize, relVerticesSize+1, relVerticesSize+2}); // Front
                    // faces.push_back({relVerticesSize+3, relVerticesSize+1, relVerticesSize+2}); // Front

                    // faces.push_back({relVerticesSize+4, relVerticesSize+5, relVerticesSize+6}); // Back
                    // faces.push_back({relVerticesSize+7, relVerticesSize+5, relVerticesSize+6}); // Back

                    // faces.push_back({relVerticesSize+1, relVerticesSize+5, relVerticesSize+3}); // Top
                    // faces.push_back({relVerticesSize+7, relVerticesSize+5, relVerticesSize+3}); // Top

                    // faces.push_back({relVerticesSize, relVerticesSize+4, relVerticesSize+6}); // Bottom
                    // faces.push_back({relVerticesSize, relVerticesSize+2, relVerticesSize+6}); // Bottom
                
                    // faces.push_back({relVerticesSize+2, relVerticesSize+6, relVerticesSize+3}); // Right
                    // faces.push_back({relVerticesSize+7, relVerticesSize+6, relVerticesSize+3}); // Right

                    // faces.push_back({relVerticesSize, relVerticesSize+4, relVerticesSize+1}); // Left
                    // faces.push_back({relVerticesSize+5, relVerticesSize+4, relVerticesSize+1}); // Left
                }
            }
        }
    }

    std::cout << vis << ", " << vertices.size() << std::endl;

    out.addVertexPositions(vertices);
    out.addFaceIndices(faces);

    out.write("output.ply", happly::DataFormat::ASCII);
}

int main(int argc, char *argv[]){
    Params p(argc, argv);

    Voxels voxels = {};
    voxels.reserve(SCENE_SIZE*SCENE_SIZE*SCENE_SIZE);
    voxels[0] = true;

    generateVoxels(voxels);
    return 0;
}