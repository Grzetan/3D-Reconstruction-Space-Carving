#include <iostream>
#include "happly.h"
#include "types.h"
#include "params.h"
#include "utils.h"
#include <map>
#include <vector>
#include "bmplib.h"

int main(int argc, char *argv[]){
    Params p(argc, argv);

    // Create voxel grid
    Voxels voxels = {};
    voxels.reserve(SCENE_SIZE*SCENE_SIZE*SCENE_SIZE);
    AABB box({0,0,0}, {VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE});

    // Read sample image
    BMP img("./images/whiteBox.bmp");

    // Camera position and rotation in real world relative to center of turn table
    Vec3 cameraPos(VOXEL_SIZE * SCENE_SIZE / 2, -1, VOXEL_SIZE * SCENE_SIZE / 2);
    Vec3 cameraDir(0, 1, 0);

    // Voxel grid center
    Vec3 gridCenter(VOXEL_SIZE * SCENE_SIZE / 2, VOXEL_SIZE * SCENE_SIZE / 2, VOXEL_SIZE * SCENE_SIZE / 2);

    // Camera parameters
    int xFOV = 70;
    int yFOV = 60;
    
    int W = img.get_width(), H = img.get_height();

    // Angles of rays shot from each pixel
    double oneTickX = (double) xFOV / (double) W;
    double oneTickY = (double) yFOV / (double) H;
    double startAngleX = (double) -xFOV / 2.0;
    double startAngleY = (double) -yFOV / 2.0;
    double angleX, angleY;

    int N_IMAGES = 1;
    double angle = 0.25*M_PI / N_IMAGES;

    Vertices v;

    rotateAroundZAxis(cameraDir, angle);

    // v.push_back({cameraDir.x, cameraDir.y, cameraDir.z});

    // rotateUsingQuaterion(cameraDir, M_PI / 4, RotationType::UP_DOWN);

    // v.push_back({cameraDir.x, cameraDir.y, cameraDir.z});

    for(int x=0; x<W; x++){
        for(int y=0; y<H; y++){
            angleX = (startAngleX + x * oneTickX);
            angleY = (startAngleY + y * oneTickY);
            Vec3 pixelDir = rotateUsingQuaterion(cameraDir, degrees2radians(angleY), RotationType::UP_DOWN);
            pixelDir = rotateUsingQuaterion(pixelDir, degrees2radians(angleX), RotationType::LEFT_RIGHT);

            v.push_back({pixelDir.x, pixelDir.y, pixelDir.z});
        }
    }

    // for(int i=0; i<N_IMAGES; i++){
    //     // Rotate camera direction and position
    //     Vec3 offset = cameraPos - gridCenter;
    //     rotateAroundZAxis(offset, angle);
    //     cameraPos = offset + gridCenter;

    //     rotateAroundZAxis(cameraDir, angle);

    //     std::cout << cameraDir.x << ", " << cameraDir.y << ", " << cameraDir.z << std::endl;
    //     std::cout << cameraPos.x << ", " << cameraPos.y << ", " << cameraPos.z << std::endl;

    //     // v.push_back({cameraDir.x, cameraDir.y, cameraDir.z});
    //     // v.push_back({cameraPos.x, cameraPos.y, cameraPos.z});

    //     for(int x=0; x<W; x++){
    //         for(int y=0; y<H; y++){
    //             Pixel pixel = img.get_pixel(x, y);
    //             // Dont remove voxels for object pixels
    //             if(pixel.r == 255 && pixel.g == 255 && pixel.b == 255) continue;

    //             angleX = (startAngleX + x * oneTickX);
    //             angleY = (startAngleY + y * oneTickY);
    //             std::cout << angleX << ", " << angleY << std::endl;
    //             // Vec3 dir{angleX + cameraDir.x, angleY + cameraDir.y, cameraDir.z};
    //             // Ray r(dir, cameraPos);

    //             // // std::cout << dir.x << ", " << dir.y << ", " << dir.z << std::endl;


    //             // rayGridTraversal(r, voxels, box);
    //         }
    //     }
    // }

    std::cout << img.get_width() << ", " << img.get_height() << std::endl;

    happly::PLYData out;
    out.addVertexPositions(v);
    out.write("test.ply");

    // generateVoxels(voxels);
    return 0;
}