#include <FreeImage.h>

#include "radiate/radiate.hh"

using Radiate::Point3f;
using Radiate::Vector3f;

static void usage(char* prgm)
{
    fprintf(stderr, "usage: %s [options]*\n"
                    "\n"
                    "Options:\n"
                    "    [-in        <input-mesh>]\n"
                    "    [-out       <output-file>]\n"
                    "    [-width     <img-pixel-width = 1920>]\n"
                    "    [-height    <img-pixel-height = 1080>]\n"
                    "    [-per-pixel <pixel-to-world = 0.001>]\n"
                    "    [-eye       <camera-origin (x, y, z) = auto]\n"
                    "    [-window    <window-origin (x, y, z) = auto]\n"
                    "    [-up        <up-vector (x, y, z) = (0, 1, 0)]\n"
                    "    [-ortho     <do-orthographic = 0>]\n",
                    prgm);
    exit(1);
}

static void parse3f(char* s, Vector3f& v)
{
    if (sscanf((const char*) s, "(%f, %f, %f)", &v[0], &v[1], &v[2]) != 3) {
        fprintf(stderr, "Unable to parse triple %s\n", s);
        exit(1);
    }
}

int main(int argc, char** argv)
{
    Radiate::SceneManager mgr;

    char* in_path = nullptr;
    char* out_path = nullptr;

    int height, width;
    mgr.getDimensions(&height, &width);

    float per_pixel;
    Point3f eye;
    Point3f window;
    Vector3f up;
    bool do_orthographic;
    mgr.getView(&per_pixel, &eye, &window, &up, &do_orthographic);

    for (int i = 1; (i + 1) < argc; i += 2) {
#define OPT(s) (strcmp(argv[i], s) == 0)
        if (OPT("-in")) {
            in_path = argv[i + 1];
        } else if (OPT("-out")) {
            out_path = argv[i + 1];
        } else if (OPT("-width")) {
            width = atoi(argv[i + 1]);
        } else if (OPT("-height")) {
            height = atoi(argv[i + 1]);
        } else if (OPT("-per-pixel")) {
            per_pixel = atof(argv[i + 1]);
        } else if (OPT("-eye")) {
            parse3f(argv[i + 1], eye);
        } else if (OPT("-window")) {
            parse3f(argv[i + 1], window);
        } else if (OPT("-up")) {
            parse3f(argv[i + 1], up);
        } else if (OPT("-ortho")) {
            do_orthographic = atoi(argv[i + 1]);
        } else {
            fprintf(stderr, "Unknown parameter: %s\n", argv[i]);
            usage(argv[0]);
        }
#undef OPT
    }

    if (!in_path || !out_path) {
        fprintf(stderr, "Missing input/output paths\n");
        usage(argv[0]);
    }

    FreeImage_Initialise(false);
    mgr.setDimensions(height, width);
    mgr.setView(per_pixel, eye, window, up, do_orthographic);
    mgr.Add(in_path);
    mgr.Render();
    mgr.getPixelBuffer()->Write(out_path);
    FreeImage_DeInitialise();

    return 0;
}
