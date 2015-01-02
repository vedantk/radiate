#pragma once

#include <string>
#include <cassert>

#include "radiate/util.hh"

namespace Radiate
{

enum kd_axes
{
    KD_X = 0,
    KD_Y = 1,
    KD_Z = 2,
    KD_NONE = 3
};

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

    KDNode() : tris(nullptr) {}

    inline bool isLeaf() { return !parent.not_leaf; }

    inline uint32_t getLeftChild() { return parent.left_index; }

    inline uint32_t getRightChild() { return parent.left_index + 1; }

    // Check if this is the last triangle in the leaf.
    static bool hasNext(Triangle* T);

    // Mark this triangle as not being the last triangle in the leaf.
    static Triangle* setNext(Triangle* T);

    // Retrieve the base triangle pointer with all tags cleared.
    static Triangle* getTriangle(Triangle* T);

    // A special destructor for leaf nodes.
    void Destroy();

#define KD_FOR_EACH(leaf, T, j)                                   \
    for (Triangle* T = KDNode::getTriangle(leaf->tris[j]);        \
         T; T = KDNode::hasNext(leaf->tris[j]) ?                  \
                KDNode::getTriangle(leaf->tris[++j]) : nullptr)
} __attribute__((packed, aligned(8)));

struct KDMetaNode
{
    KDNode kdsplit;
    BoundingSphere sphere;
};

class KDTree
{
public:
    ~KDTree();

    // Construct a KD-tree out of the given triangles.
    void Create(Mesh& tris, BoundingBox& bbox);

    // Return the first triangle the ray intersects, or nullptr.
    Triangle* Trace(Ray& ray, float* t);

private:
    uint64_t total_leaves;
    uint64_t total_leaf_tris;
    std::vector<KDMetaNode> kd_vec;

    void kd_build(uint32_t cur_idx, Mesh& tris, BoundingBox& bbox);

    Triangle* ray_trace(uint32_t cur_idx, Ray& ray, float* t);
};

}; // end namespace Radiate
