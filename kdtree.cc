#include "radiate/kdtree.hh"

#define KD_QUAD_SAMPLES     8
#define NR_SEGMENTS         (KD_QUAD_SAMPLES + 1)
#define KD_EXACT_THRESH     64

using namespace std;

namespace Radiate {

// Is there a next Triangle*?
bool KDNode::hasNext(Triangle* T)
{
    return (ptrdiff_t(T) & 2ULL) != 0;
}

// Set the "not-last-Triangle*" bit.
Triangle* KDNode::setNext(Triangle* T)
{
    return (Triangle*) (ptrdiff_t(T) | 2ULL);
}

// Remove the "not-last-Triangle*" bit.
Triangle* KDNode::getTriangle(Triangle* T)
{
    return (Triangle*) (ptrdiff_t(T) & (~2ULL));
}

// Free the Triangle** array in leaves.
void KDNode::Destroy()
{
    if (isLeaf()) {
        assert(tris);
        delete[] tris;
    }
}

KDTree::~KDTree()
{
    for (size_t i = 0; i < kd_vec.size(); ++i) {
        kd_vec[i].kdsplit.Destroy();
    }
}

void KDTree::Create(Mesh& tris)
{
    size_t nr_tris = tris.size();
    size_t sz_estimate = 4 * nr_tris;
    total_leaves = total_leaf_tris = 0;

    puts(":: Creating kd-tree...");

    kd_vec.reserve(sz_estimate);
    kd_vec.emplace_back(KDMetaNode());
    kd_build(0, tris);

    printf(":: kd-tree built with %zu nodes partitioning %zu triangles.\n",
            kd_vec.size(), nr_tris);
    printf("There are %zu leaf nodes containing a total of %zu triangles.\n",
            total_leaves, total_leaf_tris);
}

// A basic kd-tree traversal algorithm for a single ray.
Triangle* KDTree::ray_trace(uint32_t cur_idx, Ray& ray, float* t)
{
    KDMetaNode* node = &kd_vec[cur_idx];
    if (!node->sphere.Intersect(ray)) {
        return NULL;
    }

    KDNode* cur = &node->kdsplit;

    if (cur->isLeaf()) {
        Triangle* closest = NULL;
        float min_dist = INFINITY;

        size_t j = 0;
        KD_FOR_EACH(cur, T, j) {
            if (T->Intersect(ray, t) && *t < min_dist) {
                closest = T;
                min_dist = *t;
            }
        }

        *t = min_dist;
        return closest;
    }

    int split_axis = cur->parent.split_axis;
    float split_pos = cur->parent.split_pos;

    float t_boundary = (split_pos - ray.O[split_axis]) / ray.D[split_axis];

    if (ray.O[split_axis] < split_pos) {
        Triangle* L = ray_trace(cur->getLeftChild(), ray, t);

        if ((L && *t <= t_boundary) || ray.D[split_axis] <= 0) {
            return L;
        }

        return ray_trace(cur->getRightChild(), ray, t);
    } else {
        Triangle* R = ray_trace(cur->getRightChild(), ray, t);

        if ((R && *t <= t_boundary) || ray.D[split_axis] >= 0) {
            return R;
        }

        return ray_trace(cur->getLeftChild(), ray, t);
    }
}

Triangle* KDTree::Trace(Ray& ray, float* t)
{
    return ray_trace(0, ray, t);
}

// Is this Triangle to the left of the splitting plane?
static bool kd_on_left(Triangle* T, float split_pos, int axis)
{
    return T->v1[axis] < split_pos
           || T->v2[axis] < split_pos
           || T->v3[axis] < split_pos;
}

// Is this Triangle to the right of the splitting plane?
static bool kd_on_right(Triangle* T, float split_pos, int axis)
{
    return T->v1[axis] >= split_pos
           || T->v2[axis] >= split_pos
           || T->v3[axis] >= split_pos;
}

// Compute the surface-area-heuristic for the given axis split.
static float kd_eval_sah(Mesh& tris, float split, int axis)
{
    size_t lcnt = 0;
    size_t rcnt = 0;
    BoundingBox leftbb;
    BoundingBox rightbb;
    BoundingBox whole;

    for (size_t i = 0; i < tris.size(); ++i) {
        Triangle* T = tris[i];
        if (kd_on_left(T, split, axis)) {
            ++lcnt;
            leftbb.Add(T);
        }
        if (kd_on_right(T, split, axis)) {
            ++rcnt;
            rightbb.Add(T);
        }
    }

    whole.Add(leftbb);
    whole.Add(rightbb);

    float SA_l = leftbb.HalfSurfaceArea();
    float SA_r = rightbb.HalfSurfaceArea();
    float SA_w = whole.HalfSurfaceArea();
    assert((isfinite(SA_l) && SA_w >= SA_l) || !isfinite(SA_l));
    assert((isfinite(SA_r) && SA_w >= SA_r) || !isfinite(SA_r));
    assert(lcnt + rcnt >= tris.size());
    return (lcnt * SA_l/SA_w) + (rcnt * SA_r/SA_w);
}

// Minimize the SAH cost over the set of all axis splits in the mesh.
static float kd_exact_sah(Mesh& tris, float& best_split, int& best_axis,
                          BoundingBox& box)
{
    float best_sah = INFINITY;
    Vector3f widths = box.top - box.bottom;
    best_axis = max_axis(widths);

    for (size_t i = 0; i < tris.size(); ++i) {
        Triangle* T = tris[i];

        float sah = kd_eval_sah(tris, T->v1[best_axis], best_axis);
        if (sah < best_sah) {
            best_sah = sah;
            best_split = T->v1[best_axis];
            best_axis = best_axis;
        }
    }
    return best_sah;
}

// (axis-posn, SAH-Cost(axis-posn))
typedef pair<float, float> SAHSample;

static float poly_eval(float a, float b, float c, float x, float x1, float x2)
{
    return a + (x - x1)*b + (x - x1)*(x - x2)*c;
}

// Determine the SAH cost by minimizing piecewise quadratic functions.
static float kd_quad_approx(Mesh& tris, float& best_split, int& best_axis,
                            BoundingBox& box)
{
    float best_sah = INFINITY;
    for (int axis = KD_X; axis < KD_NONE; ++axis) {
        float hi = box.top[axis];
        float lo = box.bottom[axis];
        float width = hi - lo;
        float alpha = width / NR_SEGMENTS;
        float min_sah = INFINITY;
        float max_sah = -INFINITY;

        // Sample the axis interval evenly.
        vector<SAHSample> samples;
        samples.reserve(KD_QUAD_SAMPLES * 3);
        for (int i = 1; i <= KD_QUAD_SAMPLES; ++i) {
            float pos = lo + i*alpha;
            float sah = kd_eval_sah(tris, pos, axis);
            samples.push_back(make_pair(pos, sah));
            min_sah = min(min_sah, sah);
            max_sah = max(max_sah, sah);

            if (sah < best_sah) {
                best_sah = sah;
                best_split = pos;
                best_axis = axis;
            }
        }

        // Add samples wherever appreciable changes in Range{SAH} occur.
        float sah_range = max_sah - min_sah;
        float sah_thresh = sah_range / NR_SEGMENTS;
        for (int i = 0; i < NR_SEGMENTS; ++i) {
            float prev_sah = i == 0 ?
                             tris.size() : samples[i-1].second;
            float next_sah = i == KD_QUAD_SAMPLES ?
                             tris.size() : samples[i].second;
            int nr_extra = roundf(fabs(next_sah - prev_sah) / sah_thresh);

            float seg0 = lo + i*alpha;
            float beta = alpha / (nr_extra + 1);
            for (int j = 1; j <= nr_extra; ++j) {
                float pos = seg0 + j*beta;
                float sah = kd_eval_sah(tris, pos, axis);
                samples.push_back(make_pair(pos, sah));

                if (sah < best_sah) {
                    best_sah = sah;
                    best_split = pos;
                    best_axis = axis;
                }
            }
        }

        // Sort the samples by position for simpler interpolation.
        sort(samples.begin(), samples.end(),
            [](const SAHSample L, const SAHSample R)
            {
                return L.first < R.first;
            });

        // Find extrema by interpolating point triples, update the best split.
        for (size_t i = 0; i <= samples.size() - 3; ++i) {
            float x1 = samples[i].first;
            float x2 = samples[i+1].first;
            float x3 = samples[i+2].first;

            float f1 = samples[i].second;
            float f2 = samples[i+1].second;
            float f3 = samples[i+2].second;
            float a = f1;
            float b = (f2 - f1) / (x2 - x1);
            float c = (f3 - f1 - b*(x3 - x1)) / ((x3 - x1) * (x3 - x2));

            float approx_pos = 0.5 * (-b/c + x2 + x1);
            float approx_cost = poly_eval(a, b, c, approx_pos, x1, x2);

            if (x1 <= approx_pos && approx_pos <= x3 && approx_cost < best_sah)
            {
                best_sah = approx_cost;
                best_split = approx_pos;
                best_axis = axis;
            }
        }
    }
    return best_sah;
}

// Check if the mesh should be a leaf. If not, partition it.
static bool kd_mesh_split(KDNode* cur, Mesh& tris, Mesh& right_tris,
                          BoundingBox& bbox)
{
    float best_split = 0.f;
    int best_axis = KD_NONE;
    float best_sah = INFINITY;

    if (tris.size() == 1) {
        return true;
    } else if (tris.size() <= KD_EXACT_THRESH) {
        best_sah = kd_exact_sah(tris, best_split, best_axis, bbox);
    } else {
        best_sah = kd_quad_approx(tris, best_split, best_axis, bbox);
    }

    if (best_sah >= float(tris.size())) {
        return true;
    }

    cur->parent.split_axis = best_axis;
    cur->parent.split_pos = best_split;
    right_tris.reserve(tris.size() / 2);

    for (size_t i = 0; i < tris.size(); ++i) {
        Triangle* T = tris[i];
        bool on_left = kd_on_left(T, best_split, best_axis);
        bool on_right = kd_on_right(T, best_split, best_axis);

        if (on_right) {
            right_tris.push_back(T);
            if (!on_left) {
                Triangle* tail = tris.back();
                tris[i] = tail;
                tris.pop_back();
                --i;
            }
        } else {
            assert(on_left);
            continue;
        }
    }

    assert(tris.size() != 0 && right_tris.size() != 0);
    return false;
}

// Piece together internal and leaf nodes.
void KDTree::kd_build(uint32_t cur_idx, Mesh& tris)
{
    assert(tris.size() > 0);

    Mesh right_tris;
    BoundingBox bbox;
    KDMetaNode* node = &kd_vec[cur_idx];
    KDNode* cur = &node->kdsplit;

    bbox.Add(tris);
    node->sphere.Add(bbox);

    bool make_leaf = kd_mesh_split(cur, tris, right_tris, bbox);

    if (make_leaf) {
        assert(right_tris.size() == 0);

        ++total_leaves;
        total_leaf_tris += tris.size();

        cur->tris = new Triangle*[tris.size()];
        for (size_t i = 0; i < tris.size(); ++i) {
            if (i == tris.size() - 1) {
                cur->tris[i] = tris[i];
            } else {
                cur->tris[i] = KDNode::setNext(tris[i]);
            }
        }
        tris.clear();
    } else {
        assert(right_tris.size() > 0);

        cur->parent.not_leaf = 1;
        cur->parent.left_index = kd_vec.size();
        kd_vec.emplace_back(KDMetaNode());
        kd_vec.emplace_back(KDMetaNode());
        kd_build(cur->getLeftChild(), tris);
        kd_build(cur->getRightChild(), right_tris);
    }
}

}; // end namespace Radiate
