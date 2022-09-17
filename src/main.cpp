#include <iostream>
#include <filesystem>
#include <algorithm>
#include "happly.h"
#include "types.h"
#include "utils.h"
#include <map>
#include <vector>
#include "bmplib.h"
#include "argparse.hpp"

using namespace argparse;

int main(int argc, char *argv[]){
    std::cout << "Parsing arguments..." << std::endl;

    ArgumentParser args("Space Carving");

    args.add_argument("--scene_size")
           .default_value<int>(200)
           .help("Number of voxels in each dimention")
           .scan<'i', int>();

    args.add_argument("--voxel_size")
        .default_value<double>(0.2)
        .help("Voxel size in Virtual world")
        .scan<'g', double>();

    args.add_argument("--voxel_rl_size")
        .default_value<double>(1)
        .help("Voxel size in real world")
        .scan<'g', double>();

    args.add_argument("path")
        .help("Path to images");

    try {
        args.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << args;
        std::exit(1);
    }

    // Create voxel grid
    Voxels voxels(args.get<int>("--scene_size"),
                  args.get<double>("--voxel_size"),
                  args.get<double>("--voxel_rl_size"));

    AABB box({0,0,0}, {voxels.VOXEL_SIZE*voxels.SCENE_SIZE, voxels.VOXEL_SIZE*voxels.SCENE_SIZE, voxels.VOXEL_SIZE*voxels.SCENE_SIZE});

    // Camera position and rotation in real world relative to center of turn table
    Vec3 cameraPos(voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2, -voxels.SCENE_SIZE*voxels.VOXEL_SIZE, voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2);
    Vec3 cameraDir(0, 1, 0);

    // Voxel grid center
    Vec3 gridCenter(voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2, voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2, voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2);

    // Camera parameters
    double xFOV = 46.7;
    double yFOV = 46.7;
    xFOV = degrees2radians(xFOV);
    yFOV = degrees2radians(yFOV);

    std::cout << "Reading images..." << std::endl;

    // Get every image in provided path and sort them by filename
    std::string path = args.get<std::string>("path");;
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

    auto start = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    for(int i=0; i<N_IMAGES; i++){
        std::cout << "Image " << i << "/" << N_IMAGES << "\t\r" << std::flush;
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

    auto end = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    auto time_ = end.count() - start.count();

    std::cout << std::endl << "Carving took: " << time_ << ". Average time per image: " << time_ / N_IMAGES << std::endl;

    Vec3 bbox = getBoundingBox(voxels);
    std::cout << "Bounding box: (x: " << bbox.x << "mm, y: " << bbox.y << "mm, z: " << bbox.z << "mm)" << std::endl;

    std::cout << "Rendering carved space..." << std::endl;

    generateVoxels(voxels);
    return 0;
}