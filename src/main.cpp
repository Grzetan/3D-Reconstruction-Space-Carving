#include <iostream>
#include "happly.h"
#include "types.h"
#include "params.h"
#include <vector>

void generateVoxels(Voxels& voxels){
    happly::PLYData out;

    Vertices vertices;
    Faces faces;

    int x,y,z;

    for(x=0; x<SCENE_SIZE; x++){
        for(y=0; y<SCENE_SIZE; y++){
            for(z=0; z<SCENE_SIZE; z++){
                if(voxels[x + y * SCENE_SIZE + z*SCENE_SIZE*SCENE_SIZE] == false){
                    // Front, bottom, left
                    double mainX = x * VOXEL_SIZE, mainY = y * VOXEL_SIZE, mainZ = z * VOXEL_SIZE;
                    double offsetX = mainX + VOXEL_SIZE, offsetY = mainY + VOXEL_SIZE, offsetZ = mainZ + VOXEL_SIZE;

                    const size_t relVerticesSize = vertices.size();

                    // Generate vertices
                    vertices.push_back({mainX, mainY, mainZ}); // 0: Front bottom left
                    vertices.push_back({mainX, offsetY, mainZ}); // 1: Front top left
                    vertices.push_back({offsetX, mainY, mainZ}); // 2: Front bottom right
                    vertices.push_back({offsetX, offsetY, mainZ}); // 3: Front top right
                    vertices.push_back({mainX, mainY, offsetZ}); // 4: Back bottom left
                    vertices.push_back({mainX, offsetY, offsetZ}); // 5: Back top left
                    vertices.push_back({offsetX, mainY, offsetZ}); // 6: Back bottom right
                    vertices.push_back({offsetX, offsetY, offsetZ}); // 7: Back top right

                    // Generate faces
                    faces.push_back({relVerticesSize, relVerticesSize+1, relVerticesSize+2}); // Front
                    faces.push_back({relVerticesSize+3, relVerticesSize+1, relVerticesSize+2}); // Front

                    faces.push_back({relVerticesSize+4, relVerticesSize+5, relVerticesSize+6}); // Back
                    faces.push_back({relVerticesSize+7, relVerticesSize+5, relVerticesSize+6}); // Back

                    faces.push_back({relVerticesSize+1, relVerticesSize+5, relVerticesSize+3}); // Top
                    faces.push_back({relVerticesSize+7, relVerticesSize+5, relVerticesSize+3}); // Top

                    faces.push_back({relVerticesSize, relVerticesSize+4, relVerticesSize+6}); // Bottom
                    faces.push_back({relVerticesSize, relVerticesSize+2, relVerticesSize+6}); // Bottom
                
                    faces.push_back({relVerticesSize+2, relVerticesSize+6, relVerticesSize+3}); // Right
                    faces.push_back({relVerticesSize+7, relVerticesSize+6, relVerticesSize+3}); // Right

                    faces.push_back({relVerticesSize, relVerticesSize+4, relVerticesSize+1}); // Left
                    faces.push_back({relVerticesSize+5, relVerticesSize+4, relVerticesSize+1}); // Left
                }
            }
        }
    }

    out.addVertexPositions(vertices);
    out.addFaceIndices(faces);

    out.write("output.ply", happly::DataFormat::Binary);
}

int main(int argc, char *argv[]){
    Params p(argc, argv);

    Voxels voxels = {};
    voxels.reserve(SCENE_SIZE*SCENE_SIZE*SCENE_SIZE);
    voxels[0] = true;

    generateVoxels(voxels);
    return 0;
}