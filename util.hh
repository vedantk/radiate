#pragma once

#include <cmath>
#include <string>
#include <limits>
#include <cstdio>
#include <cstdint>

#include "radiate/aux.hh"

namespace Radiate
{

#define ZEROV Point3f(0.0, 0.0, 0.0)
#define INFTYV Point3f(INFINITY, INFINITY, INFINITY)
#define EPSILON 1e-8

void xerr(std::string s);

void printv(Vector3f& v);

void printvln(Vector3f& v);

static inline float square(float x)
{
    return x * x;
}

static inline float dsquare(Vector3f& v)
{
    return v.dot(v);
}

static inline int pick_median(int i, int j)
{
    return i + ((j - i) / 2);
}

static inline float fequal(float x, float y)
{
    return fabs(x - y) < EPSILON;
}

static inline int max_axis(Vector3f& v)
{
    if (v[0] > v[1]) {
        return (v[0] > v[2]) ? 0 : 2;
    } else {
        return (v[1] > v[2]) ? 1 : 2;
    }
}

struct Ray
{
    Point3f O;
    Vector3f D; // ||D|| = 1

    Ray(Point3f origin, Vector3f direction)
        : O(origin),
          D(direction.normalized())
    {}

    inline Point3f getPoint(float t) { return O + t*D; }
};

struct Triangle
{
    Point3f v1, v2, v3;
    Vector3f e1; // v2 - v1
    Vector3f e2; // v3 - v1
    Vector3f normal; // e1 x e2, s.t ||n|| = 1

    // Initialize the given triangle.
    static void Init(Triangle* T, Point3f _v1, Point3f _v2, Point3f _v3);

    // Allocate a triangle from an object pool, then initialize it.
    static Triangle* Create(Point3f _v1, Point3f _v2, Point3f _v3);

    bool Intersect(Ray& ray, float* t);
};

struct BoundingBox
{
    Point3f top;
    Point3f bottom;

    BoundingBox()
        : top(-INFTYV),
          bottom(INFTYV)
    {}

    void Add(Point3f& pt);

    void Add(Triangle* T);

    void Add(Mesh& mesh);

    void Add(BoundingBox& box);

    // Compute δxδy + δxδz + δyδz.
    float HalfSurfaceArea();
};

void printb(BoundingBox& bbox);

void printbln(BoundingBox& bbox);

struct BoundingSphere
{
    Point3f center;
    float rsquared;

    BoundingSphere()
        : center(ZEROV),
          rsquared(0.0)
    {}

    void Add(BoundingBox& bbox);

    bool Intersect(Ray& ray);
};

}; // end namespace Radiate
