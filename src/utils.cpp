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
 
    if(tmin < 0 || tmax < 0) return false;

    tstart = tmin;
    tend = tmax;

    return true; 
}

// https://www.youtube.com/watch?v=lJdEX3w0xaY
// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.3443&rep=rep1&type=pdf
void rayGridTraversal(Ray& ray, Voxels& voxels, const AABB& box){
    double tstart, tend;

    // If ray doesn't hit AABB of voxels, skip
    if(!rayAABBIntersection(ray, box, tstart, tend)){
        return;
    }

    // Calculate enter point and voxel
    Vec3 enterPoint = ray.orig + ray.dir * tstart;

    // Get binary direction of ray
    int xDir = ray.dir.x > 0 ? 1 : ray.dir.x < 0 ? -1 : 0;
    int yDir = ray.dir.y > 0 ? 1 : ray.dir.y < 0 ? -1 : 0;
    int zDir = ray.dir.z > 0 ? 1 : ray.dir.z < 0 ? -1 : 0;

    int enterVoxel[3];

    enterVoxel[0] = std::floor(enterPoint.x / (VOXEL_SIZE + -xDir * 1e-6));
    enterVoxel[1] = std::floor(enterPoint.y / (VOXEL_SIZE + -yDir * 1e-6));
    enterVoxel[2] = std::floor(enterPoint.z / (VOXEL_SIZE + -zDir * 1e-6));

    // Get world position of grid boundaries ray is going to hit
    double xBound = (xDir > 0 ? enterVoxel[0]+1 : enterVoxel[0]) * VOXEL_SIZE;
    double yBound = (yDir > 0 ? enterVoxel[1]+1 : enterVoxel[1]) * VOXEL_SIZE;
    double zBound = (zDir > 0 ? enterVoxel[2]+1 : enterVoxel[2]) * VOXEL_SIZE;

    // Get the exact time of intersection with each boundary
    double xt = std::abs((xBound - enterPoint.x) / (ray.dir.x + 1e-6));
    double yt = std::abs((yBound - enterPoint.y) / (ray.dir.y + 1e-6));
    double zt = std::abs((zBound - enterPoint.z) / (ray.dir.z + 1e-6));

    // Time needed to travel through voxel for each axis
    double xDelta = VOXEL_SIZE * xDir / (ray.dir.x + 1e-6);
    double yDelta = VOXEL_SIZE * yDir / (ray.dir.y + 1e-6);
    double zDelta = VOXEL_SIZE * zDir / (ray.dir.z + 1e-6);

    // What is the termination index for each axis
    int xOut = xDir < 0 ? -1 : SCENE_SIZE;
    int yOut = yDir < 0 ? -1 : SCENE_SIZE;
    int zOut = zDir < 0 ? -1 : SCENE_SIZE;

    xt = (xt <= 0) ? 100000 : xt;
    yt = (yt <= 0) ? 100000 : yt;
    zt = (zt <= 0) ? 100000 : zt;

    // Traverse until exit
    while(true){
        voxels[idx(enterVoxel[0], enterVoxel[1], enterVoxel[2])] = true;

        if(xt < yt && xt < zt){
            enterVoxel[0] += xDir;
            if(enterVoxel[0] == xOut) break;
            xt += xDelta;
        }else if(yt < zt){
            enterVoxel[1] += yDir;
            if(enterVoxel[1] == yOut) break;
            yt += yDelta;
        }else{
            enterVoxel[2] += zDir;
            if(enterVoxel[2] == zOut) break;
            zt += zDelta;
        }
    }
}

void rotateAroundZAxis(Vec3& v, double angle){
    double l = v.len();
    double x = v.x * std::cos(angle) - v.y * std::sin(angle);
    double y = v.x * std::sin(angle) + v.y * std::cos(angle);
    v.x = x;
    v.y = y;
    v.normalize();
    v = v * l;
}

double degrees2radians(double degrees){
    return degrees * (M_PI / 180.0);
}

Vec3 rotateUsingQuaterion(Vec3& v, double angle, RotationType type){
    // Find axis to rotate about (Vec3 is initialized with 0s)
    Vec3 axis;
    if(type == RotationType::LEFT_RIGHT){
        axis.z = 1;
    }else if(type == RotationType::UP_DOWN){
        axis.x = v.y;
        axis.y = -v.x;
    }

    Quaternion R(angle, axis);
    Quaternion P(v);
    Quaternion invR = R.inverse();

    Quaternion newP = R.hamilton(P).hamilton(invR);

    return {newP.b, newP.c, newP.d};
}
