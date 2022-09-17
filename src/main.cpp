#include <iostream>
#include <filesystem>
#include <algorithm>
#include "happly.h"
#include "types.h"
#include "params.h"
#include "utils.h"
#include <map>
#include <vector>
#include "bmplib.h"

int main(int argc, char *argv[]){
    std::cout << "Parsing params ..." << std::endl;

    Params p(argc, argv);

    // Create voxel grid
    Voxels voxels = {};
    voxels.reserve(SCENE_SIZE*SCENE_SIZE*SCENE_SIZE);
    AABB box({0,0,0}, {VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE, VOXEL_SIZE*SCENE_SIZE});

    // Camera position and rotation in real world relative to center of turn table
    Vec3 cameraPos(VOXEL_SIZE * SCENE_SIZE / 2, -SCENE_SIZE*VOXEL_SIZE, VOXEL_SIZE * SCENE_SIZE / 2);
    Vec3 cameraDir(0, 1, 0);

    // Voxel grid center
    Vec3 gridCenter(VOXEL_SIZE * SCENE_SIZE / 2, VOXEL_SIZE * SCENE_SIZE / 2, VOXEL_SIZE * SCENE_SIZE / 2);

    // Camera parameters
    double xFOV = 46.7;
    double yFOV = 46.7;
    xFOV = degrees2radians(xFOV);
    yFOV = degrees2radians(yFOV);

    std::cout << "Reading images..." << std::endl;

    // Get every image in provided path and sort them by filename
    std::string path = "./out";
    auto dir = std::filesystem::directory_iterator(path);
    std::vector<std::string> images = {};

    for (const auto & entry : dir)
        images.push_back(entry.path());

    struct Comp
    {
        inline int splitByDot(std::string& s){
            std::istringstream iss(s);
            std::string token;
            std::getline(iss, token, '.');

            return stoi(token);
        }

        inline bool operator() (const std::string& s1, const std::string& s2)
        {
            std::string filename1 = std::filesystem::path(s1).filename();
            std::string filename2 = std::filesystem::path(s2).filename();

            return splitByDot(filename1) < splitByDot(filename2);
        }
    };

    std::sort(images.begin(), images.end(), Comp());

    // Read sample image to get width and height
    BMP sampleImg(images[0]);

    int W = sampleImg.get_width(), H = sampleImg.get_height();

    // Angles of rays shot from each pixel
    double oneTickX = xFOV / (double) W;
    double oneTickY = yFOV / (double) H;
    double startAngleX = -xFOV / 2.0;
    double startAngleY = -yFOV / 2.0;
    double angleX, angleY;

    int N_IMAGES = images.size();
    double angle = 2*M_PI / N_IMAGES;

    std::cout << "Carving space..." << std::endl;

    for(int i=0; i<N_IMAGES; i++){
        std::cout << "Image " << i << std::endl;
        // Read current image
        BMP img(images[i]);

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
                if(pixel.r != 0 && pixel.g != 0 && pixel.b != 0) continue;

                angleX = (startAngleX + x * oneTickX);
                angleY = (startAngleY + y * oneTickY);
                Vec3 pixelDir = rotateUsingQuaterion(cameraDir, angleY, RotationType::UP_DOWN);
                pixelDir = rotateUsingQuaterion(pixelDir, angleX, RotationType::LEFT_RIGHT);
                Ray r(pixelDir, cameraPos);

                rayGridTraversal(r, voxels, box);
            }
        }
    }

    std::cout << "Rendering carved space" << std::endl;

    generateVoxels(voxels);
    return 0;
}