#include "utils.h"

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

size_t vIdx(Key key, UniqueVertices& uniqueVertices){
    if(uniqueVertices.find(key) == uniqueVertices.end()){
        uniqueVertices[key] = uniqueVertices.size();
    }

    return uniqueVertices[key];
}

void generateVoxels(Voxels& voxels){
    happly::PLYData out;

    Vertices vertices;
    Faces faces;
    UniqueVertices uniqueVertices;

    int x,y,z;

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
                        faceVisibility[CubeFaces(i)] = isCubeBackground(x+cx, y+cy, z+cz, voxels);
                    }

                    // Front, bottom, left
                    double mainX = x * VOXEL_SIZE, mainY = y * VOXEL_SIZE, mainZ = z * VOXEL_SIZE;
                    double offsetX = mainX + VOXEL_SIZE, offsetY = mainY + VOXEL_SIZE, offsetZ = mainZ + VOXEL_SIZE;

                    // Generate only visible faces
                    if(faceVisibility[LEFT]){
                        Key keyA = {mainX, offsetY, mainZ};
                        Key keyB = {mainX, mainY, mainZ};
                        Key keyC = {mainX, mainY, offsetZ};
                        Key keyD = {mainX, offsetY, offsetZ};

                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyD, uniqueVertices), vIdx(keyC, uniqueVertices)});
                    }

                    if(faceVisibility[RIGHT]){
                        Key keyA = {offsetX, offsetY, mainZ};
                        Key keyB = {offsetX, mainY, mainZ};
                        Key keyC = {offsetX, mainY, offsetZ};
                        Key keyD = {offsetX, offsetY, offsetZ};

                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyD, uniqueVertices), vIdx(keyC, uniqueVertices)});
                    }

                    if(faceVisibility[TOP]){
                        Key keyA = {mainX, offsetY, mainZ};
                        Key keyB = {offsetX, offsetY, mainZ};
                        Key keyC = {mainX, offsetY, offsetZ};
                        Key keyD = {offsetX, offsetY, offsetZ};

                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                        faces.push_back({vIdx(keyD, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                    }

                    if(faceVisibility[BOTTOM]){
                        Key keyA = {mainX, mainY, mainZ};
                        Key keyB = {offsetX, mainY, mainZ};
                        Key keyC = {mainX, mainY, offsetZ};
                        Key keyD = {offsetX, mainY, offsetZ};

                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                        faces.push_back({vIdx(keyD, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                    }

                    if(faceVisibility[FRONT]){
                        Key keyA = {mainX, mainY, mainZ};
                        Key keyB = {mainX, offsetY, mainZ};
                        Key keyC = {offsetX, mainY, mainZ};
                        Key keyD = {offsetX, offsetY, mainZ};

                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                        faces.push_back({vIdx(keyD, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                    }

                    if(faceVisibility[BACK]){
                        Key keyA = {mainX, mainY, offsetZ};
                        Key keyB = {mainX, offsetY, offsetZ};
                        Key keyC = {offsetX, mainY, offsetZ};
                        Key keyD = {offsetX, offsetY, offsetZ};

                        faces.push_back({vIdx(keyA, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                        faces.push_back({vIdx(keyD, uniqueVertices), vIdx(keyB, uniqueVertices), vIdx(keyC, uniqueVertices)});
                    }
                }
            }
        }
    }

    double x1, y1, z1;
    vertices.reserve(uniqueVertices.size());
    for(int i=0; i<uniqueVertices.size(); i++){
        vertices.push_back({0,0,0});
    }

    UniqueVertices::iterator it;
    for(it=uniqueVertices.begin(); it!=uniqueVertices.end(); it++){
        std::tie(x1,y1,z1) = it->first;
        vertices[it->second] = {x1, y1, z1};
    }

    out.addVertexPositions(vertices);
    out.addFaceIndices(faces);

    out.write("output.ply", happly::DataFormat::ASCII);
}

bool rayAABBIntersection(const Ray& ray,
                         const AABB& box,
                         double& tstart,
                         double& tend) {
    float tmin, tmax, tymin, tymax, tzmin, tzmax; 
 
    tmin = (box.bounds[ray.sign[0]].x - ray.orig.x) * ray.invDir.x; 
    tmax = (box.bounds[1-ray.sign[0]].x - ray.orig.x) * ray.invDir.x; 
    tymin = (box.bounds[ray.sign[1]].y - ray.orig.y) * ray.invDir.y; 
    tymax = (box.bounds[1-ray.sign[1]].y - ray.orig.y) * ray.invDir.y; 
 
    if ((tmin > tymax) || (tymin > tmax)) 
        return false; 
 
    if (tymin > tmin) 
        tmin = tymin; 
    if (tymax < tmax) 
        tmax = tymax; 
 
    tzmin = (box.bounds[ray.sign[2]].z - ray.orig.z) * ray.invDir.z; 
    tzmax = (box.bounds[1-ray.sign[2]].z - ray.orig.z) * ray.invDir.z; 
 
    if ((tmin > tzmax) || (tzmin > tmax)) 
        return false; 
 
    if (tzmin > tmin) 
        tmin = tzmin; 
    if (tzmax < tmax) 
        tmax = tzmax; 
 
    tstart = tmin;
    tend = tmax;

    return true; 
}