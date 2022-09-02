#include "params.h"

Params::Params(int argc, char *argv[]){
    // TODO: Handle exceptions and missing values
    for(int i=0; i<argc; i++){
        if(argv[i] == "-cubeSize"){
            voxelCubeSize = atoi(argv[i+1]);
        }
    }
}