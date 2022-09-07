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
    Vec3 cameraPos(VOXEL_SIZE * SCENE_SIZE / 2, -20, VOXEL_SIZE * SCENE_SIZE / 2);
    Vec3 cameraDir(0, 1, 0);

    // Voxel grid center
    Vec3 gridCenter(VOXEL_SIZE * SCENE_SIZE / 2, VOXEL_SIZE * SCENE_SIZE / 2, VOXEL_SIZE * SCENE_SIZE / 2);

    // Camera parameters
    double xFOV = 55;
    double yFOV = 55;
    xFOV = degrees2radians(xFOV);
    yFOV = degrees2radians(yFOV);

    int W = img.get_width(), H = img.get_height();

    // Angles of rays shot from each pixel
    double oneTickX = xFOV / (double) W;
    double oneTickY = yFOV / (double) H;
    double startAngleX = -xFOV / 2.0;
    double startAngleY = -yFOV / 2.0;
    double angleX, angleY;

    int N_IMAGES = 36;
    double angle = 2*M_PI / N_IMAGES;

    for(int i=0; i<N_IMAGES; i++){
        // Rotate camera direction and position
        Vec3 offset = cameraPos - gridCenter;
        rotateAroundZAxis(offset, angle);
        cameraPos = offset + gridCenter;

        rotateAroundZAxis(cameraDir, angle);

        // For every pixel calculate it's direction and travere grid
        for(int x=0; x<W; x++){
            for(int y=0; y<H; y++){
                Pixel pixel = img.get_pixel(x, y);
                // Dont remove voxels for object pixels
                if(pixel.r == 255 && pixel.g == 255 && pixel.b == 255) continue;

                angleX = (startAngleX + x * oneTickX);
                angleY = (startAngleY + y * oneTickY);
                Vec3 pixelDir = rotateUsingQuaterion(cameraDir, angleY, RotationType::UP_DOWN);
                pixelDir = rotateUsingQuaterion(pixelDir, angleX, RotationType::LEFT_RIGHT);
                Ray r(pixelDir, cameraPos);

                rayGridTraversal(r, voxels, box);
            }
        }
    }

    generateVoxels(voxels);
    return 0;
}