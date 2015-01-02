#include <algorithm>

#include <assimp/scene.h>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>

#include <FreeImage.h>

#include "radiate/kdtree.hh"
#include "radiate/radiate.hh"

namespace Radiate
{

PixelBuffer::~PixelBuffer()
{
    FIBITMAP* fip = (FIBITMAP*) buf;
    if (fip) {
        FreeImage_Unload(fip);
    }
}

void PixelBuffer::Resize(int h, int w)
{
    FIBITMAP* fip = (FIBITMAP*) buf;
    if (fip) {
        FreeImage_Unload(fip);
    }

    fip = FreeImage_Allocate(w, h, 32, FI_RGBA_RED_MASK, FI_RGBA_BLUE_MASK,
            FI_RGBA_GREEN_MASK);
    assert(FreeImage_GetLine(fip) / FreeImage_GetWidth(fip) == sizeof(pixel));
    if (!fip) {
        xerr("Could not allocate bitmap");
    }

    buf = (void*) fip;
}

struct pixel* PixelBuffer::getScanline(int y)
{
    FIBITMAP* fip = (FIBITMAP*) buf;
    return (struct pixel*) FreeImage_GetScanLine(fip, y);
}

void PixelBuffer::Write(char* out_path)
{
    FIBITMAP* fip = (FIBITMAP*) buf;
    if (!FreeImage_Save(FIF_PNG, fip, (const char*) out_path, 0)) {
        xerr("Failed to save image");
    }
}

SceneManager::SceneManager()
    : width(1920),
      height(1080),
      focal(0.05),
      per_pixel(0.0),
      eye(ZEROV),
      window(ZEROV),
      up(Vector3f(0, 1, 0)),
      do_orthographic(false),
      kdtree(nullptr),
      bbox(nullptr),
      in(ZEROV),
      side(ZEROV)
{
    bbox = new BoundingBox;
    pixbuf.Resize(height, width);
}

SceneManager::~SceneManager()
{
    if (kdtree) {
        delete kdtree;
    }
    if (bbox) {
        delete bbox;
    }
}

// Convert an aiFace vertex into a Point3f.
static Point3f toPoint3f(aiMesh* mesh, aiFace* face, int vtx)
{
    aiVector3D v = mesh->mVertices[face->mIndices[vtx]];
    return Point3f(v.x, v.y, v.z);
}

// Extract all triangles from the model.
static void fill_mesh(Mesh& tris, const aiScene* scene, BoundingBox* bbox)
{
    if (!scene->mNumMeshes) {
        xerr("No objects found in mesh");
    }

    for (size_t i = 0; i < scene->mNumMeshes; ++i) {
        aiMesh* mesh = scene->mMeshes[i];

        for (size_t j = 0; j < mesh->mNumFaces; ++j) {
            aiFace* face = &mesh->mFaces[j];

            if (face->mNumIndices != 3) {
                xerr("Found a non-triangular face -- unsupported");
            }

            auto T = Triangle::Create(toPoint3f(mesh, face, 0),
                                      toPoint3f(mesh, face, 1),
                                      toPoint3f(mesh, face, 2));
            tris.push_back(T);
            bbox->Add(T);
        }
    }
}

void SceneManager::AddMesh(char* in_path)
{
    // Load the scene, then deallocate it.
    {
        printf(":: Loading mesh %s...\n", in_path);
        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile((const char*) in_path,
                                                 aiProcess_Triangulate);
        if (!scene) {
            xerr("Error loading scene -- importer.ReadFile");
        }

        mesh.reserve(mesh.size() + scene->mNumMeshes);
        fill_mesh(mesh, scene, bbox);
    }
    
    // Set the viewing eye and window if a good one doesn't exist.
    if (!checkViewConsistency()) {
        float Bx = bbox->bottom[0];
        float By = bbox->bottom[1];
        float Bz = bbox->bottom[2];
        float Wx = bbox->top[0] - bbox->bottom[0];
        float Wy = bbox->top[1] - bbox->bottom[1];

        window[0] = eye[0] = Bx + (Wx / 2.0); // Midpoint{bbox.X}
        window[1] = eye[1] = By + (Wy / 2.0); // Midpoint{bbox.Y}

        // Let the desired angle-of-view be φ. To find Ez, set:
        // tan(φ) = (Wy/2)/(Bz-Ez) => Ez = Bz-(Wy/2)/tan(φ), with φ = π/5.
        float phi = M_PI / 5.0;
        eye[2] = Bz - ((Wy / 2.0) / tan(phi));
        window[2] = eye[2] + focal * (Bz - eye[2]);

        // Now find the size of the viewing window. Start with φ:
        // tan(φ) = (dh/2)/(Wz - Ez) => dh = 2(Wz - Ez)tan(φ)
        per_pixel = 2.0*(window[2] - eye[2])*tan(phi)/height;
    }

    if (!checkViewConsistency()) {
        xerr("unable to correct viewing frustrum");
    }
}

void SceneManager::getDimensions(int* h, int* w)
{
    *w = width;
    *h = height;
}

void SceneManager::setDimensions(int h, int w)
{
    if (h < 0 || w < 0) {
        xerr("Image dimensions must be positive");
    }

    width = w;
    height = h;
    pixbuf.Resize(h, w);
}

bool SceneManager::checkViewConsistency()
{
    in = window - eye;
    side = up.cross(in);

    return in.norm() > 0.0
        && fequal(in.dot(up), 0.0)
        && fequal(up.dot(side), 0.0)
        && fequal(side.dot(in), 0.0)
        && !fequal(focal, 0.0)
        && per_pixel > 0.0;
}

void SceneManager::getView(float* f, float* pp, Point3f* _eye,
        Point3f* _window, Vector3f* _up, bool* _do_ortho)
{
    *f = focal;
    *pp = per_pixel;
    *_eye = eye;
    *_window = window;
    *_up = up;
    *_do_ortho = do_orthographic;
}

void SceneManager::setView(float f, float pp, Point3f& _eye,
        Point3f& _window, Vector3f& _up, bool _do_ortho)
{
    focal = f;
    per_pixel = pp;
    eye = _eye;
    window = _window;
    up = _up;
    do_orthographic = _do_ortho;
}

void SceneManager::getBoundingBox(Point3f* bottom, Point3f* top)
{
    *bottom = bbox->bottom;
    *top = bbox->top;
}

void SceneManager::Render()
{
    if (!mesh.size()) {
        xerr("No objects to render");
    }

    if (!checkViewConsistency()) {
        xerr("inconsistent viewing frustrum");
    }

    if (!kdtree) {
        kdtree = new KDTree;
        kdtree->Create(mesh, *bbox);
        assert(!mesh.size());
    }

    in.normalize();
    side.normalize();
    up.normalize();

    float vert_size = height * per_pixel;
    float horiz_size = width * per_pixel;

    Vector3f half_vert = (vert_size / 2) * up;
    Vector3f half_horiz = (horiz_size / 2) * side;
    Point3f TL = window + half_vert - half_horiz;
    Vector3f horiz_step = (horiz_size / width) * side;
    Vector3f vert_step = (-vert_size / height) * up;

    printf("--------------------------------------------------------------\n");
    printf("Scene configuration:\n");
    printf(" - in     = "); printvln(in);
    printf(" - side   = "); printvln(side);
    printf(" - up     = "); printvln(up);
    printf(" - eye    = "); printvln(eye);
    printf(" - window = "); printvln(window);
    printf(" - bbox   = "); printbln(*bbox);
    printf(" - win.TL = "); printvln(TL);
    printf(" - win.up = "); printvln(half_vert);
    printf(" - win.-> = "); printvln(half_horiz);
    printf(" - f/frac = %.9f\n", focal);
    printf(" - d/pix  = %.9f\n", per_pixel);
    printf("--------------------------------------------------------------\n");

    float last_pctg = 0.0;
    for (int r = 0; r < height; ++r) {
        float pctg = (100.0 * r) / height;
        if (pctg - last_pctg >= 1.0) {
            printf("Render: %.3f%% done.\n", pctg);
            last_pctg = pctg;
        }

        struct pixel* line = pixbuf.getScanline(height - r - 1);

#pragma omp parallel for
        for (int c = 0; c < width; ++c) {
            Point3f cell = TL + (c * horiz_step) + (r * vert_step);
            Point3f& source = do_orthographic ? window : cell;
            Ray ray = Ray(cell, source - eye);

            float t;
            Triangle* T = kdtree->Trace(ray, &t);
            if (T) {
                Vector4f alight(0, 0, 0, 255);
                Vector4f dlight(0, 0, 0, 255);
                Vector4f slight(0, 100, 0, 255);
                Vector3f Lm = Vector3f(0, 0, -1);
                float T_Lm = T->normal.dot(Lm);
                float ambient = 0.25;
                float diffuse = 0.55 * T_Lm;
                Vector3f Rm = 2.0 * T_Lm * T->normal - Lm;
                float specular = powf(Rm.dot(-ray.D), 2.0);
                Vector4f color = (ambient * alight) + (diffuse * dlight)
                               + (specular * slight);

                line[c].r = uint8_t(std::min(255.f, color[0]));
                line[c].g = uint8_t(std::min(255.f, color[1]));
                line[c].b = uint8_t(std::min(255.f, color[2]));
                line[c].a = 255;
            } else {
                line[c].r = line[c].g = line[c].b = 0;
                line[c].a = 255;
            }
        }
    }
}

PixelBuffer* SceneManager::getPixelBuffer()
{
    return &pixbuf;
}

}; // end namespace Radiate
