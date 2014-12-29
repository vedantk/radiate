#pragma once

#include "radiate/aux.hh"

namespace Radiate
{

struct pixel
{
    uint8_t r, g, b, a;
};

class KDTree;
class BoundingBox;

class PixelBuffer
{
public:
    PixelBuffer()
        : buf(nullptr)
    {}

    ~PixelBuffer();

    // Resize the buffer.
    void Resize(int h, int w);

    // Get the r-th row from the bottom of the image.
    struct pixel* getRow(int r);

    // Write the image.
    void Write(char* out_path);

private:
    void* buf;
};

class SceneManager
{
public:
    SceneManager();

    ~SceneManager();

    // Add a new mesh into the scene.
    void Add(char* in_path);

    // Get the dimensions (in pixels) of the final image.
    void getDimensions(int* h, int* w);

    // Set the dimensions (in pixels) of the final image.
    void setDimensions(int h, int w);

    // Check if the viewing frustrum is consistent.
    bool checkViewConsistency();

    // Find the current viewing frustrum configuration.
    void getView(float* pp, Point3f* _eye, Point3f* _window, Vector3f* _up,
            bool* _do_ortho);

    // Optional viewing frustrum configuration. This overrides any
    // automatically selected viewing parameters.
    void setView(float pp, Point3f& _eye, Point3f& _window, Vector3f& _up,
            bool _do_ortho);

    // Find the bounding box of the scene, s.t top-bottom is the box diagonal.
    void getBoundingBox(Point3f* bottom, Point3f* top);

    // Render the scene.
    void Render();

    // Get a (h x w)-sized image representing the finished scene.
    PixelBuffer* getPixelBuffer();

private:
    int width;
    int height;
    float per_pixel;
    Point3f eye;
    Point3f window;
    Vector3f up;
    bool do_orthographic;
    PixelBuffer pixbuf;
    KDTree* kdtree; 
    Mesh mesh;
    BoundingBox* bbox;
};

}; // end namespace Radiate
