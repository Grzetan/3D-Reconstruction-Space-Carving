
# 3D reconstruction using Space Carving

Implementation of space (voxel) carving written in C++. 

---

# Table of contents

1. [Requirements](#requirements)
2. [How to use](#how-to-use)
3. [How it works](#how-it-works)
4. [Paramaters](#parameters)
    - [Path](#path)
    - [`--scene_size`](#scene_size)
    - [`--voxel_size`](#voxel_size)
    - [`--voxel_rl_size`](#voxel_rl_size)
    - [`--segmentation_thresh`](#thresh)
5. [Examples](#examples)
6. [Tutorial on generating test datasets in blender](#blender)


---

# Requirements <a name="requirements"></a>

- Linux (tested on Ubuntu 20.04 LTS)

- g++ compiler
```bash
sudo apt intall g++
```

- CMake
```bash
sudo snap install cmake
```

---

# How to use <a name="how-to-use"></a>

Go into main directory and type
```bash
cmake . && make
```
After that new binary called `main` will be created. Call it using
```bash
./main ./path/to/images
```

Please read `Parameters` section down below to use program correctly.


---

# How it works <a name="how-it-works"></a>

Space Carving is pretty simple algorithm for 3D reconstruction. I recommend those two videos to grasp the main idea behind this algorithm.
https://www.youtube.com/watch?v=XpjWXDY8aes
https://www.youtube.com/watch?v=cGs90KF4oTc
The program takes in grayscale `.bmp` images of one object taken from different angles. Something like this:
![img](./gallery/8.png)
and turns it into something like this:
![img](./gallery/8.png)
It does so by first segmenting an image and shooting rays for every pixel and removing those voxels which get hit by `background ray`.
`background ray` is a ray that corresponds to pixel qualified as background during segmentation phase.

Program contains 3 main steps:

### Image segmentation
First each image in input directory will have to be segmented. It means that each pixel can be qualified as either background or object. Once it segments an image, it know which rays will be carving the space (background pixels) and which will be skipped (object pixels).

### Space Carving
In this phase program will shoot rays from different angles that will piece by piece carve the space and hopefully photographed object will be created.

### Generating model
The last step is generating a `.ply` file with created object. Only visible faces will be generated to make this object as light as possible.

---

# Parameters <a name="parameters"></a>

## Path (Mandatory) <a name="path"></a>

Path to folder with `.bmp` images. Number of images in folder can be random. Please note that names should be called `1.bmp`, `2.bmp`, `3.bmp` ..., otherwise program won't execute correctly.
```bash
./main ./path/to/images
```
---

## `--scene_size` (Not required) <a name="scene_size"></a>

Specify the size of starting voxel grid. Remember that the bigger this argument is the slower and more precise program will be. I find value `200` to be good trade off between quality and execution time.

Default value is 200 so starting voxel grid would be `200x200x200`.

```bash
./main path/to/images --scene_size 300
```

## `--voxel_size` (Not required) <a name="voxel_size"></a>

Size of voxel in generated model. (It does not have any impact on algorithm).

Default value is `0.2`. If you want your object to be bigger or smaller, adjust this value.

```bash
./main path/to/images --voxel_size 0.5
```
---

## `--voxel_rl_size` (Not required) <a name="voxel_rl_size"></a>

Size of one voxel in milimeters. To measure it you can space carve simple object (like a cube) of known size and than divide the real life length (or any other dimention) by number of voxels that create this dimention. 

Default value is `1`

```bash
./main path/to/images --voxel_rl_size 4
```
---

## `--segmentation_thresh` (Not required) <a name="thresh"></a>

Threshold for segmenting an input image. If pixel is brighter than `segmentation_thresh` it will be qualified as part of an object.

Default value is `60`. Range: `<1, 255>`.

Example:
```bash
./main path/to/images --segmentation_thresh 100
```

# Gallery <a name="gallery"></a>


# Generate datasets in blender <a name="blender"></a>

TO create dataset with blender:

1. Set camera pos using 'G' and entering desired translation

2. Set camera rotation using 'R' and entering desired rotation

3. Change background color by clicking red globe on the right side

4. Set up output frames, and format (BMP) using printer icon on the right side

5. Make object emit light using red sphere icon on the right sight and setting emit color to white.

6. After loading object, click 'N' in "Rotation secting" below Z axis, enter #frame/X where X controls rotation speed (5.73 is 10degrees each frame).

7. Generate dataset using 'Render' on the top left corner.
