#include "api.h"

OutCylinder spaceCarve(const char* path, unsigned int sceneSize = 100, double rlVoxelSize = 1, double xFOV = 46.7, double yFOV = 46.7, int thresh = 60){
       // Create voxel grid
    Voxels voxels(sceneSize, 0.2, 1);

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

    // Get every image in provided path and sort them by filename
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

    if(N_IMAGES % 4 == 0){
        std::cout << "Adjusting rotation axis..." << std::endl;
        int degrees90Rotation = N_IMAGES / 4;
        std::vector<std::array<int, 2>> objectHorizontalSize = {}; // Array for storing size of detected object in X axis

        // For 4 images (0, 90, 180 and 270) calculate the size of image. This will be necessary to calculate rotation axis
        for(int i=0; i<4; i++){
            BMP img(images[i*degrees90Rotation]);
            segmentImage(img);
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
    segmentImage(sampleImg);
    calibrateVoxels(sampleImg, voxels, cameraPos, gridCenter, oneTickX, startAngleX);

    const int THRESH = thresh;

    // if(args.get<bool>("--save_segmented")){
    //     std::filesystem::remove_all("segmented_images");
    //     std::filesystem::create_directory("segmented_images");
    // }

    for(int i=0; i<N_IMAGES; i++){
        std::cout << "Image " << i << "/" << N_IMAGES << "\t\r" << std::flush;
        // Read current image
        BMP img(images[i]);

        segmentImage(img);

        // if(args.get<bool>("--save_segmented")){
        //     img.write("segmented_images/output" + std::to_string(i) + ".bmp");
        // }
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

    // if(args.get<bool>("--filter")){
    //     std::cout << "Removing small groups of voxels..." << std::endl;
    //     removeSingleVoxels(voxels);
    // }

    // Remove gate
    // Vector of areas to remove (x1,y1,x2,y2)
    std::vector<std::array<int, 6>> areasToRemove = { {50,80,60, 85,120,200}, {70,70,160, 200,200,190}, {130,80,60, 200,120,200}, {60,90,60, 150,120,100} };
    removeGate(voxels, areasToRemove);

    Bbox bbox = getBoundingBox(voxels);

    Cylinder cylinder = getCylinder(voxels, bbox, axisXPos, axisZPos);

    return {cylinder.r, cylinder.h};
}