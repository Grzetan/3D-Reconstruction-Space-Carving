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
    Vec3 cameraPos(VOXEL_SIZE * SCENE_SIZE / 2 - 0.01, 0, VOXEL_SIZE * SCENE_SIZE / 2 - 0.01);
    Vec3 cameraDir(0, 1, 0);

    // Camera parameters
    int xFOV = 70;
    int yFOV = 40;
    
    // For every background pixel (black) shoot ray at corresponding angle
    int W = img.get_width(), H = img.get_height();

    double oneTickX = (double) xFOV / (double) W;
    double oneTickY = (double) yFOV / (double) H;
    double startAngleX = (double) -xFOV / 2.0;
    double startAngleY = (double) -yFOV / 2.0;
    double angleX, angleY;

    for(int x=0; x<W; x++){
        for(int y=0; y<H; y++){
            Pixel pixel = img.get_pixel(x, y);
            // Dont remove voxels for object pixels
            if(pixel.r == 255 && pixel.g == 255 && pixel.b == 255) continue;

            angleX = (startAngleX + x * oneTickX) / 90.0;
            angleY = (startAngleY + y * oneTickY) / 90.0;
            Vec3 dir{angleX + cameraDir.x, cameraDir.y, angleY + cameraDir.z};
            Ray r(dir, cameraPos);

            rayGridTraversal(r, voxels, box);
        }
    }

    std::cout << img.get_width() << ", " << img.get_height() << std::endl;

    generateVoxels(voxels);
    return 0;
}