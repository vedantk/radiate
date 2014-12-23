#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>

#include <boost/pool/poolfwd.hpp>
#include <boost/pool/object_pool.hpp>

#include "math.hh"

static void usage(char* prgm)
{
    fprintf(stderr, "usage: %s [options]*\n"
                    "    [-in       <input-mesh>]\n"
                    "    [-out      <output-file = out.png>]\n"
                    "    [-width    <img-pixel-width = 1920>]\n"
                    "    [-height   <img-pixel-height = 1080>]\n"
                    "    [-camera   <camera-origin (x, y, z) = auto]\n"
                    "    [-window   <window-origin (x, y, z) = auto]\n"
                    );
    exit(1);
}

int main(int argc, char** argv)
{
    
}
