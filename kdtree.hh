#pragma once

#include <string>
#include <cassert>

#include "util.hh"

namespace Radiate
{

struct Ray
{
    Point3f O;
    Vector3f D; // ||D|| = 1

    Ray(Point3f origin, Vector3f direction)
        : O(origin), D(direction.normalized())
    {}

    Point3f getPoint(float t) { return O + t*D; }
};

struct Triangle
{
    Point3f v1, v2, v3;
    Vector3f e1; // v2 - v1
    Vector3f e2; // v3 - v1
    Vector3f normal; // e1 x e2

    static Triangle* Create(Point3f _v1, Point3f _v2, Point3f _v3);
    static void Init(Triangle* T, Point3f _v1, Point3f _v2, Point3f _v3);

    bool Intersect(Ray& ray, float* t);
};

enum kd_axes
{
    KD_X = 0,
    KD_Y = 1,
    KD_Z = 2,
    KD_NONE = 3
};

typedef std::vector<Triangle*> Mesh;

union KDNode
{
    // Parent nodes: store the splitting plane and left child index.
    struct
    {
        float split_pos;
        uint32_t split_axis : 2;
        uint32_t left_index : 29;
        uint32_t not_leaf : 1;
    } parent;

    // Leaf nodes: store an array of Triangle pointers with link tags.
    Triangle** tris;

    KDNode() : tris(NULL) {}

    bool isLeaf() { return !parent.not_leaf; }
    uint32_t getLeftChild() { return parent.left_index; }
    uint32_t getRightChild() { return parent.left_index + 1; }

    static bool hasNext(Triangle* T);
    static Triangle* setNext(Triangle* T);
    static Triangle* getTriangle(Triangle* T);

    void Destroy();

#define KD_FOR_EACH(leaf, T, j)                               \
    for (Triangle* T = KDNode::getTriangle(leaf->tris[j]);    \
         T; T = KDNode::hasNext(leaf->tris[j]) ?              \
                KDNode::getTriangle(leaf->tris[++j]) : NULL)
};

struct BoundingSphere
{
    Point3f C;
    float rsquared;

    // Check if a ray enters the bounding sphere.
    bool Intersect(Ray& ray);
};

struct KDMetaNode
{
    KDNode node;
    BoundingSphere sphere;
};

class KDTree
{
public:
    ~KDTree();

    void Create(Mesh& tris);

    // Shoot a ray into the scene. Return a triangle if one is hit.
    Triangle* Trace(Ray& ray, float* t);

private:
    uint64_t total_leaves;
    uint64_t total_leaf_tris;
    std::vector<KDMetaNode> kd_vec;

    void kd_build(uint32_t cur_idx, Mesh& tris);
    Triangle* ray_trace(uint32_t cur_idx, Ray& ray, float* t);
};

}; // end namespace Radiate
