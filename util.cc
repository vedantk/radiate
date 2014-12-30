#include "radiate/util.hh"

#include <boost/pool/poolfwd.hpp>
#include <boost/pool/object_pool.hpp>

namespace Radiate
{

using namespace std;

void xerr(std::string s)
{
    fprintf(stderr, "fatal: %s\n", s.c_str());
    exit(1);
}

void printv(Vector3f& v)
{
    printf("(%.3f, %.3f, %.3f)", v[0], v[1], v[2]);
}

void printvln(Vector3f& v)
{
    printv(v);
    printf("\n");
}

static boost::object_pool<Triangle> TrianglePool;

void Triangle::Init(Triangle* T, Point3f _v1, Point3f _v2, Point3f _v3)
{
    T->v1 = _v1;
    T->v2 = _v2;
    T->v3 = _v3;
    T->e1 = _v2 - _v1;
    T->e2 = _v3 - _v1;
    T->normal = T->e1.cross(T->e2).normalized();
}

// Check that the triangle vertices are finite numbers.
static bool normal_triangle(Triangle* tri)
{
    for (int i = 0; i < 3; ++i) {
        if (!isfinite(tri->v1[i])) return false;
        if (!isfinite(tri->v2[i])) return false;
        if (!isfinite(tri->v3[i])) return false;
    }
    return true;
}

Triangle* Triangle::Create(Point3f _v1, Point3f _v2, Point3f _v3)
{
    Triangle* T = TrianglePool.malloc();
    Triangle::Init(T, _v1, _v2, _v3);

    if (normal_triangle(T)) {
        return T;
    } else {
        xerr("Triangle::Create -- abnormal triangle detected");
        return nullptr;
    }
}

// A standard Möller-Trumbore ray-triangle intersection test.
bool Triangle::Intersect(Ray& ray, float* t)
{
    Vector3f P, Q, T;
    float det, inv_det, u, v;

    P = ray.D.cross(e2);
    det = e1.dot(P);
    if (det > -EPSILON && det < EPSILON) {
        return false;
    }
    inv_det = 1.0 / det;

    T = ray.O - v1;

    u = T.dot(P) * inv_det;
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    Q = T.cross(e1);

    v = ray.D.dot(Q) * inv_det;
    if (v < 0.0 || (u + v > 1.0)) {
        return false;
    }

    *t = e2.dot(Q) * inv_det;

    return *t > EPSILON;
}

void BoundingBox::Add(Point3f& pt)
{
    top[0] = max(top[0], pt[0]);
    top[1] = max(top[1], pt[1]);
    top[2] = max(top[2], pt[2]);
    bottom[0] = min(bottom[0], pt[0]);
    bottom[1] = min(bottom[1], pt[1]);
    bottom[2] = min(bottom[2], pt[2]);
}

void BoundingBox::Add(Triangle* T)
{
    // Using .cwise{Max,Min} makes a bottleneck, so unroll it:
    top[0] = max(max(top[0], T->v1[0]), max(T->v2[0], T->v3[0]));
    top[1] = max(max(top[1], T->v1[1]), max(T->v2[1], T->v3[1]));
    top[2] = max(max(top[2], T->v1[2]), max(T->v2[2], T->v3[2]));
    bottom[0] = min(min(bottom[0], T->v1[0]), min(T->v2[0], T->v3[0]));
    bottom[1] = min(min(bottom[1], T->v1[1]), min(T->v2[1], T->v3[1]));
    bottom[2] = min(min(bottom[2], T->v1[2]), min(T->v2[2], T->v3[2]));
}

void BoundingBox::Add(Mesh& mesh)
{
    for (size_t i = 0; i < mesh.size(); ++i) {
        Add(mesh[i]);
    }
}

void BoundingBox::Add(BoundingBox& box)
{
    Add(box.top);
    Add(box.bottom);
}

float BoundingBox::HalfSurfaceArea()
{
    Point3f ext = top - bottom;
    return ext[0]*ext[1] + ext[0]*ext[2] + ext[1]*ext[2];
}

void BoundingSphere::Add(BoundingBox& bbox)
{
    Vector3f H = 0.5 * (bbox.top - bbox.bottom);
    center = bbox.bottom + H;
    rsquared = dsquare(H);
}

bool BoundingSphere::Intersect(Ray& ray)
{
    // <O + tD - C, O + tD - C> = r^2
    // t^2D.D - 2tD.C + O.O + 2O.tD - 2O.C - C.C - r^2 = 0
    Point3f d = ray.O - center;
    float b = 2.0 * d.dot(ray.D);
    float c = d.dot(d) - rsquared;
    float discriminant = square(b) - 4*c; // ||D|| = 1 => t^2 coeff = 1.
    if (discriminant < 0) {
        return false;
    }

    float e = sqrtf(discriminant);
    float t0 = (-b + e) / 2.0;
    if (t0 >= 0) {
        return true;
    }

    float t1 = (-b - e) / 2.0;
    return t1 >= 0;
}

}; // end namespace Radiate
