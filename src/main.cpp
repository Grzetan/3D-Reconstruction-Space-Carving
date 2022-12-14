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

    args.add_argument("--segmentation_thresh")
        .default_value<int>(50)
        .help("Color of an object")
        .scan<'i', int>();

    args.add_argument("--marching_cubes")
           .help("Use marching cubes algorithm when generating output PLY")
           .default_value(false)
           .implicit_value(true);

    args.add_argument("--filter")
        .help("Delete small groups of voxels that could be not hit be removing ray")
        .default_value(false)
        .implicit_value(true);

    args.add_argument("--adjust_rotation")
        .help("Automaticly adjust rotation axis based on images. Number of images in folder must be divisible by 4.")
        .default_value(false)
        .implicit_value(true);

    args.add_argument("--save_segmented")
        .help("Save segmented images to output directory")
        .default_value(false)
        .implicit_value(true);

    args.add_argument("--remove_gate")
        .help("Remove gate from output image")
        .default_value(false)
        .implicit_value(true);

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

    constexpr double CAMERA_DIST = 0.9;

    // Camera position and rotation in real world relative to center of turn table
    Vec3 cameraPos(voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2, (-CAMERA_DIST)*voxels.SCENE_SIZE*voxels.VOXEL_SIZE, voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2);
    Vec3 cameraDir(0, 1, 0);

    // Voxel grid center
    Vec3 gridCenter(voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2, voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2, voxels.VOXEL_SIZE * voxels.SCENE_SIZE / 2);

    // Camera parameters
    double xFOV = 55.7;
    double yFOV = 46.7;
    // double xFOV = 60;
    // double yFOV = 60;
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

    // If possible, adjust rotation axis
    int N_IMAGES = images.size();

    int VERTICAL_MARGIN = 200; // Margin on the top and bottom (This allows to adjust axis even if rotary table is visible)

    double axisXPos, axisZPos;

    // Variables that are passed to segment image to calcualte object brightness
    double intensitySum = 0;
    unsigned long pixelNum = 0;

    if(N_IMAGES % 4 == 0 && args.get<bool>("--adjust_rotation")){
        std::cout << "Adjusting rotation axis..." << std::endl;
        int degrees90Rotation = N_IMAGES / 4;
        std::vector<std::array<int, 2>> objectHorizontalSize = {}; // Array for storing size of detected object in X axis

        // For 4 images (0, 90, 180 and 270) calculate the size of image. This will be necessary to calculate rotation axis
        for(int i=0; i<4; i++){
            BMP img(images[i*degrees90Rotation]);
            segmentImage(img, intensitySum, pixelNum);
            int min = W, max = 0;
            for(int x=0; x<W; x++){
                for(int y=VERTICAL_MARGIN; y<H; y++){
                    Pixel pixel = img.get_pixel(x, y);
                    if(pixel.r == 0 && pixel.g == 0 && pixel.b == 0){
                        min = std::min(min, x);
                        max = std::max(max, x);
                    }
                }
            }
            objectHorizontalSize.push_back({min, max});
        }

        // Adjust position of rotation axis in X
        axisXPos = (double)(objectHorizontalSize[0][0] + (objectHorizontalSize[2][1] - objectHorizontalSize[0][0]) / 2) / W;
        axisZPos = (double)(objectHorizontalSize[1][0] + (objectHorizontalSize[3][1] - objectHorizontalSize[1][0]) / 2) / W;

        std::cout << axisXPos << ", " << axisZPos << std::endl;

        cameraPos.x = voxels.VOXEL_SIZE * voxels.SCENE_SIZE * axisXPos;
        cameraPos.z = voxels.VOXEL_SIZE * voxels.SCENE_SIZE * axisZPos;
    }

    // Angles of rays shot from each pixel
    double oneTickX = xFOV / (double) W;
    double oneTickY = yFOV / (double) H;
    double startAngleX = -xFOV / 2.0;
    double startAngleY = -yFOV / 2.0;
    double angleX, angleY;

    double angle = 2*M_PI / N_IMAGES;

    
    // Calculate RL voxel size from first image
    segmentImage(sampleImg, intensitySum, pixelNum);
    calibrateVoxels(sampleImg, voxels, cameraPos, gridCenter, oneTickX, startAngleX);

    std::cout << "Carving space..." << std::endl;

    auto start = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    const int THRESH = args.get<int>("--segmentation_thresh");

    if(args.get<bool>("--save_segmented")){
        std::filesystem::remove_all("segmented_images");
        std::filesystem::create_directory("segmented_images");
    }

    intensitySum = 0;
    pixelNum = 0;

    for(int i=0; i<N_IMAGES; i++){
        std::cout << "Image " << i << "/" << N_IMAGES << "\t\r" << std::flush;
        // Read current image
        BMP img(images[i]);

        segmentImage(img, intensitySum, pixelNum);

        if(args.get<bool>("--save_segmented")){
            img.write("segmented_images/output" + std::to_string(i) + ".bmp");
        }
        // Rotate camera direction and position
        if(i){
            Vec3 offset = cameraPos - gridCenter;
            rotateAroundZAxis(offset, angle);
            cameraPos = offset + gridCenter;

            rotateAroundZAxis(cameraDir, angle);
        }


        // For every pixel calculate it's direction and travere grid
        for(int x=0; x<W; x++){
            for(int y=0; y<H; y++){
                Pixel pixel = img.get_pixel(x, y);
                // Dont remove voxels for object pixels
                if(pixel.r == 0 && pixel.g == 0 && pixel.b == 0) continue;

                angleX = (startAngleX + x * oneTickX);
                angleY = (startAngleY + y * oneTickY);
                Vec3 pixelDir = rotateUsingQuaterion(cameraDir, angleY, RotationType::UP_DOWN);
                pixelDir = rotateUsingQuaterion(pixelDir, angleX, RotationType::LEFT_RIGHT);   
                Ray r(pixelDir, cameraPos);

                rayGridTraversal(r, voxels, box);
            }
        }
    }

    if(args.get<bool>("--filter")){
        std::cout << "Removing small groups of voxels..." << std::endl;
        removeSingleVoxels(voxels);
    }

    auto end = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );

    auto time_ = end.count() - start.count();

    std::cout << std::endl << "Carving took: " << time_ << ". Average time per image: " << time_ / N_IMAGES << std::endl;

    // Remove gate
    if(args.get<bool>("--remove_gate")){
        // Vector of areas to remove (x1,y1,x2,y2)
        std::vector<std::array<int, 6>> areasToRemove = { {50,80,60, 85,120,200}, {70,70,160, 200,200,190}, {130,80,60, 200,120,200}, {60,90,60, 150,120,100} };
        removeGate(voxels, areasToRemove);
    }

    Bbox bbox = getBoundingBox(voxels);
    std::cout << "Bounding box: (width: " << (bbox.max.x - bbox.min.x) * voxels.VOXEL_RL_SIZE << "cm, height: " << (bbox.max.z - bbox.min.z) * voxels.VOXEL_RL_SIZE << "cm, depth: " << (bbox.max.y - bbox.min.y) * voxels.VOXEL_RL_SIZE << "cm)" << std::endl;

    Cylinder cylinder = getCylinder(voxels, bbox, axisXPos, axisZPos);
    std::cout << "Cylinder: radius = " << cylinder.r << "cm, height = " << cylinder.h << "cm, center = (" << cylinder.center.x << ", " << cylinder.center.y << ", " << cylinder.center.z << ")" << std::endl;

    double objectBrightness = intensitySum / pixelNum;
    std::cout << "Object brightness: " << objectBrightness << " / 255" << std::endl;

    std::cout << "Rendering carved space..." << std::endl;

    if(args.get<bool>("--marching_cubes")){
        generateVoxelsMarchingCubes(voxels);
    }else{
        generateVoxels(voxels);
    }
    
    return 0;
}