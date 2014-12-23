#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <Eigen/Dense>

using Eigen::Vector3f;
using Eigen::Vector4f;

typedef Vector3f Point3f;

#define ZEROV Point3f(0.f, 0.f, 0.f)
#define INFTYV Point3f(INFINITY, INFINITY, INFINITY)
#define EPSILON 1e-8

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
