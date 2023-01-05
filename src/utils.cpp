#include "utils.h"

int idx(int x, int y, int z, Voxels& voxels){
    x = std::min(voxels.SCENE_SIZE-1, std::max(0, x));
    y = std::min(voxels.SCENE_SIZE-1, std::max(0, y));
    z = std::min(voxels.SCENE_SIZE-1, std::max(0, z));

    return x + y*voxels.SCENE_SIZE + z*voxels.SCENE_SIZE*voxels.SCENE_SIZE;
}

bool getVoxelValue(Voxels& voxels, int x, int y, int z){
    if(x < 0 || x >= voxels.SCENE_SIZE || y < 0 || y >= voxels.SCENE_SIZE || z < 0 || z >= voxels.SCENE_SIZE){
        return true;
    }

    return voxels.data[x + y*voxels.SCENE_SIZE + z*voxels.SCENE_SIZE*voxels.SCENE_SIZE];
}

bool isCubeBackground(int x, int y, int z, Voxels& voxels){
    if(x<0 || x>voxels.SCENE_SIZE-1 || y<0 || y>voxels.SCENE_SIZE-1 || z<0 || z>voxels.SCENE_SIZE-1){
        return true;
    }

    return voxels.data[idx(x, y, z, voxels)];
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

    for(x=0; x<voxels.SCENE_SIZE; x++){
        for(y=0; y<voxels.SCENE_SIZE; y++){
            for(z=0; z<voxels.SCENE_SIZE; z++){
                // Draw voxel if val is true and there is at least one neighbour not visible
                if(voxels.data[idx(x, y, z, voxels)] == false){
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
                    double mainX = x * voxels.VOXEL_SIZE, mainY = y * voxels.VOXEL_SIZE, mainZ = z * voxels.VOXEL_SIZE;
                    double offsetX = mainX + voxels.VOXEL_SIZE, offsetY = mainY + voxels.VOXEL_SIZE, offsetZ = mainZ + voxels.VOXEL_SIZE;

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

void generateVoxelsMarchingCubes(Voxels& voxels){
    happly::PLYData out;

    Vertices vertices;
    Faces faces;
    UniqueVertices uniqueVertices;

    int x,y,z;

    for(x=0; x<voxels.SCENE_SIZE; x++){
        for(y=0; y<voxels.SCENE_SIZE; y++){
            for(z=0; z<voxels.SCENE_SIZE; z++){
                if(voxels.data[idx(x, y, z, voxels)] == false){
                    int shapeIdx = 0;

                    // Find combination idx
                    for(int i=0; i<GUIDE.size(); i++){
                        int cx, cy, cz;
                        std::tie(cx, cy, cz) = GUIDE[i];
                        if(isCubeBackground(x+cx, y+cy, z+cz, voxels)){
                            for(int activeVertex : ACTIVE_VERTICES[i]){
                                shapeIdx |= 1 << activeVertex;
                            }
                        }
                    }

                    if(shapeIdx == 0 || shapeIdx == 255) continue;


                    // vertices.push_back({x*voxels.VOXEL_SIZE,y*voxels.VOXEL_SIZE,z*voxels.VOXEL_SIZE});
                    // vertices.push_back({x*voxels.VOXEL_SIZE+voxels.VOXEL_SIZE, y*voxels.VOXEL_SIZE, z*voxels.VOXEL_SIZE});
                    // vertices.push_back({x*voxels.VOXEL_SIZE, y*voxels.VOXEL_SIZE+voxels.VOXEL_SIZE, z*voxels.VOXEL_SIZE});
                    // faces.push_back({vertices.size() - 3, vertices.size() - 2, vertices.size() - 1});
                    for(int i=0; i<5; i++){
                        if(triTable[shapeIdx][i*3] == -1) break;
                        for(int j=0; j<3; j++){
                            vertices.push_back({x*voxels.VOXEL_SIZE+edgeTable[triTable[shapeIdx][i*3+j]][0]*voxels.VOXEL_SIZE, 
                                                y*voxels.VOXEL_SIZE+edgeTable[triTable[shapeIdx][i*3+j]][1]*voxels.VOXEL_SIZE, 
                                                z*voxels.VOXEL_SIZE+edgeTable[triTable[shapeIdx][i*3+j]][2]*voxels.VOXEL_SIZE});
                        }
                        faces.push_back({vertices.size() - 3, vertices.size() - 2, vertices.size() - 1});
                    }
                }
            }
        }
    }

    // double x1, y1, z1;
    // vertices.reserve(uniqueVertices.size());
    // for(int i=0; i<uniqueVertices.size(); i++){
    //     vertices.push_back({0,0,0});
    // }

    // UniqueVertices::iterator it;
    // for(it=uniqueVertices.begin(); it!=uniqueVertices.end(); it++){
    //     std::tie(x1,y1,z1) = it->first;
    //     vertices[it->second] = {x1, y1, z1};
    // }

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

    enterVoxel[0] = std::floor(enterPoint.x / (voxels.VOXEL_SIZE + -xDir * 1e-6));
    enterVoxel[1] = std::floor(enterPoint.y / (voxels.VOXEL_SIZE + -yDir * 1e-6));
    enterVoxel[2] = std::floor(enterPoint.z / (voxels.VOXEL_SIZE + -zDir * 1e-6));

    // Get world position of grid boundaries ray is going to hit
    double xBound = (xDir > 0 ? enterVoxel[0]+1 : enterVoxel[0]) * voxels.VOXEL_SIZE;
    double yBound = (yDir > 0 ? enterVoxel[1]+1 : enterVoxel[1]) * voxels.VOXEL_SIZE;
    double zBound = (zDir > 0 ? enterVoxel[2]+1 : enterVoxel[2]) * voxels.VOXEL_SIZE;

    // Get the exact time of intersection with each boundary
    double xt = std::abs((xBound - enterPoint.x) / (ray.dir.x + 1e-6));
    double yt = std::abs((yBound - enterPoint.y) / (ray.dir.y + 1e-6));
    double zt = std::abs((zBound - enterPoint.z) / (ray.dir.z + 1e-6));

    // Time needed to travel through voxel for each axis
    double xDelta = voxels.VOXEL_SIZE * xDir / (ray.dir.x + 1e-6);
    double yDelta = voxels.VOXEL_SIZE * yDir / (ray.dir.y + 1e-6);
    double zDelta = voxels.VOXEL_SIZE * zDir / (ray.dir.z + 1e-6);

    // What is the termination index for each axis
    int xOut = xDir < 0 ? -1 : voxels.SCENE_SIZE;
    int yOut = yDir < 0 ? -1 : voxels.SCENE_SIZE;
    int zOut = zDir < 0 ? -1 : voxels.SCENE_SIZE;

    xt = (xt <= 0) ? 100000 : xt;
    yt = (yt <= 0) ? 100000 : yt;
    zt = (zt <= 0) ? 100000 : zt;

    // Traverse until exit
    while(true){
        voxels.data[idx(enterVoxel[0], enterVoxel[1], enterVoxel[2], voxels)] = true;

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

Bbox getBoundingBox(Voxels& voxels){
    int x,y,z;
    double minX=voxels.SCENE_SIZE, maxX=0, minY=voxels.SCENE_SIZE, maxY=0, minZ=voxels.SCENE_SIZE, maxZ=0;

    for(x=0; x<voxels.SCENE_SIZE; x++){
        for(y=0; y<voxels.SCENE_SIZE; y++){
            for(z=0; z<voxels.SCENE_SIZE; z++){
                if(voxels.data[idx(x,y,z,voxels)]) continue;

                if(x < minX) minX = x;
                if(x > maxX) maxX = x;

                if(y < minY) minY = y;
                if(y > maxY) maxY = y;

                if(z < minZ) minZ = z;
                if(z > maxZ) maxZ = z;
            }
        }
    }

    return {
        {minX, minY, minZ},
        {maxX+1, maxY+1, maxZ+1}
    };
}

Cylinder getCylinder(Voxels& voxels, Bbox& bbox, const double axisXPos, const double axisZPos){    
    double axisX = (axisXPos * voxels.SCENE_SIZE);
    double axisZ = (axisZPos * voxels.SCENE_SIZE);

    std::cout << axisX << ", " << axisZ << std::endl;

    // axisX = axisZ = 100;

    // Find the longest distance between voxel and roation axis

    // int x,y,z;
    // for(x=0; x<voxels.SCENE_SIZE; x++){
    //     for(y=0; y<voxels.SCENE_SIZE; y++){
    //         for(z=0; z<voxels.SCENE_SIZE; z++){
    //             if(voxels.data[idx(x,y,z,voxels)]) continue;

    //             double distX = (x * voxels.VOXEL_SIZE - axisX);
    //             double distY = (y * voxels.VOXEL_SIZE - axisZ);
    //             double dist = sqrt(distX*distX + distY*distY);
    //             if(dist > maxDist) maxDist = dist;
    //         }
    //     }
    // }


    // Return the distance between the farthest bbox point and rotation axis

    double W = bbox.getWidth();
    double D = bbox.getDepth();
    const std::vector<std::array<double, 2>> points = {{bbox.min.x, bbox.min.y + D/2}, {bbox.min.x + W/2, bbox.min.y}, {bbox.min.x + W/2, bbox.max.y}, {bbox.max.x, bbox.min.y + D/2}};
    double maxDist = 0;

    for(auto& p : points){
        double dist = sqrt(std::pow(p[0]-axisX, 2) + std::pow(p[1]-axisZ, 2));
        if(maxDist < dist) maxDist = dist;
    }

    double r = maxDist * voxels.VOXEL_RL_SIZE;

    // Very simple way to do this

    // double w = (bbox.max.x - bbox.min.x);
    // double d = (bbox.max.y - bbox.min.y);
    // double r = (w > d) ? w / 2 : d / 2;

    return {r, (bbox.max.z - bbox.min.z) * voxels.VOXEL_RL_SIZE, {0 / voxels.VOXEL_SIZE * voxels.VOXEL_RL_SIZE, bbox.min.y, 0 / voxels.VOXEL_SIZE * voxels.VOXEL_RL_SIZE}};
}

// Create voxel group
void createGroupVoxels(Voxels& voxels, int x, int y, int z, std::vector<VoxelArea>& groups){
    VoxelArea area{{x,y,z}, {x+1, y+1, z+1}};

    int i,j;
    bool expandedX = true, expandedY = true, expandedZ = true;

    while(expandedX || expandedY || expandedZ){
        expandedX = false;
        expandedY = false;
        expandedZ = false;

        // Check if there are any voxels on "end" back perimiter (z const)
        for(i=area.start.x; i<=area.end.x; i++){
            for(j=area.start.y; j<=area.start.y; j++){
                if(!getVoxelValue(voxels, i,j,area.end.z)){
                    area.end.z++;
                    expandedZ = true;
                }
            }
        }

        // Check if there are any voxels on "end" right perimiter (x const)
        for(i=area.start.y; i<=area.end.y; i++){
            for(j=area.start.z; j<=area.start.z; j++){
                if(!getVoxelValue(voxels, area.end.x,i,j)){
                    area.end.x++;
                    expandedX = true;
                }
            }
        }

        // Check if there are any voxels on "end" bottom perimiter (y const)
        for(i=area.start.x; i<=area.end.x; i++){
            for(j=area.start.z; j<=area.start.z; j++){
                if(!getVoxelValue(voxels, i,area.end.y,j)){
                    area.end.y++;
                    expandedY = true;
                }
            }
        }
    }

    groups.push_back(area);
}

void createGroupImage(BMP& image, int x, int y, std::vector<ImageArea>& groups){
    ImageArea area{{x,y}, {x+1, y+1}};

    int i;
    bool expandedX = true, expandedY = true;

    while(expandedX || expandedY){
        expandedX = false;
        expandedY = false;

        // Check if there are any object pixels to the right
        for(i=area.start[1]; i<=area.end[1]; i++){
            if(area.end[0] >= image.get_width() || i >= image.get_height()) break;

            Pixel p = image.get_pixel(area.end[0], i);
            if(!isPixelBackground(p)){
                area.end[0]++;
                expandedX = true;
            }
        }

        // Check if there are any object pixels at the bottom
        for(i=area.start[0]; i<=area.end[0]; i++){
            if(i >= image.get_width() || area.end[1] >= image.get_height()) break;

            Pixel p = image.get_pixel(i,area.end[1]);
            if(!isPixelBackground(p)){
                area.end[1]++;
                expandedY = true;
            }
        }
    }

    groups.push_back(area);
}

// Remove voxels small groups of voxels so they don't ruin "out cylinder"
void removeSingleVoxels(Voxels& voxels){
    int x,y,z;

    std::vector<VoxelArea> detectedGroups = {};

    for(x=0; x<voxels.SCENE_SIZE; x++){
        for(y=0; y<voxels.SCENE_SIZE; y++){
            for(z=0; z<voxels.SCENE_SIZE; z++){
                if(voxels.data[idx(x,y,z,voxels)]) continue;

                // If voxel is already in group, skip
                bool isInGroup = false;
                for(auto& group : detectedGroups){
                    if(x >= group.start.x && y >= group.start.y && z >= group.start.z && 
                       x <= group.end.x && y <= group.end.y && z <= group.end.z){
                        if(!group.isValid()) voxels.data[idx(x,y,z,voxels)] = true;
                        isInGroup = true;
                        break;
                    }
                }
                if(isInGroup) continue;
                

                createGroupVoxels(voxels, x,y,z, detectedGroups);
                if(!detectedGroups[detectedGroups.size() - 1].isValid()) voxels.data[idx(x,y,z,voxels)] = true;
            }
        }
    }
}

// WHITE = background
// BLACK = object
void segmentImage(BMP& image, double& intensitySum, unsigned long& pixelNum){
    int x, y;
    std::vector<ImageArea> detectedGroups = {};

    for(x=0; x<image.get_width(); x++){
        for(y=0; y<image.get_height(); y++){
            Pixel pixel = image.get_pixel(x, y);
            // If background, change color and skip
            if(isPixelBackground(pixel) || x > 520 || y < 100){
                image.set_pixel(x, y, 255, 255, 255, 255);
                continue;
            }else{
                intensitySum += (pixel.r + pixel.g + pixel.b) / 3;
                pixelNum++;
            }
            // image.set_pixel(x, y, 0,0,0, 255);
            // continue;
            // If pixel is already in group, skip
            bool isInGroup = false;
            for(auto& group : detectedGroups){
                if(x >= group.start[0] && y >= group.start[1] && x <= group.end[0] && y <= group.end[1]){
                    // If group is not valid, change color to background (white), if valid change to object (black)
                    if(!group.isValid()){
                        image.set_pixel(x, y, 255, 255, 255, 255);
                    }else{
                        image.set_pixel(x, y, 0, 0, 0, 255);
                    }
                    isInGroup = true;
                    break;
                }
            }
            if(isInGroup) continue;

            createGroupImage(image, x, y, detectedGroups);
            if(!detectedGroups[detectedGroups.size() - 1].isValid()) image.set_pixel(x, y, 255, 255, 255, 255);
        }
    }
}

void removeGate(Voxels& voxels, std::vector<std::array<int, 6>>& areasToRemove){    
    for(auto& area : areasToRemove){
        // Remove every voxel in area
        for(int x=area[0]; x<area[3]; x++){
            for(int y=area[1]; y<area[4]; y++){
                for(int z=area[2]; z<area[5]; z++){
                    voxels.data[idx(x,y,z,voxels)] = true;
                }
            }
        }
    }
}

void calibrateVoxels(BMP& img, Voxels& voxels, Vec3& cameraPos, Vec3& gridCenter, double oneTick, double startAngle){
    int yOffset = 390;
    int min=img.get_width(), max=0;

    for(int x=0; x<img.get_width(); x++){
        Pixel pixel = img.get_pixel(x, yOffset);
        if(pixel.r == 0 && pixel.g == 0 && pixel.b == 0){
            min = (x < min) ? x : min;
            max = (x > max) ? x : max;
        }
    }
    if(min == max) std::runtime_error("Cannot detect frame");

    double angleMin = (startAngle + min * oneTick);
    double angleMax = (startAngle + max * oneTick);

    Vec3 minDir = {0,1,0};
    Vec3 maxDir = {0,1,0};

    rotateAroundZAxis(minDir, angleMin);
    rotateAroundZAxis(maxDir, angleMax);

    // std::cout << minDir.x << ", " << minDir.y << ", " << minDir.z << std::endl;
    // std::cout << maxDir.x << ", " << maxDir.y << ", " << maxDir.z << std::endl;

    minDir = minDir * cameraPos.dist(gridCenter);
    maxDir = maxDir * cameraPos.dist(gridCenter);
    double distBetween = minDir.dist(maxDir);
    double voxelsInside = distBetween / voxels.VOXEL_SIZE;
    double voxelRLSize = 28.0 / voxelsInside;
    voxels.VOXEL_RL_SIZE = voxelRLSize;
}

bool isPixelBackground(Pixel& p){
    return (p.g > p.r && p.g > p.b && p.g > 100 && p.g - p.r > 30 && p.g - p.b > 30) || (p.r < 30 && p.g < 30 && p.b < 30);
    // return (p.r < 30 && p.g < 30 && p.b < 30);
}