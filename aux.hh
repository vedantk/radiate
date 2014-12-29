#pragma once

#include <vector>

#include <eigen3/Eigen/Dense>

namespace Radiate
{

using Eigen::Vector3f;
using Eigen::Vector4f;

typedef Vector3f Point3f;

class Triangle;

typedef std::vector<Triangle*> Mesh;

}; // end namespace Radiate
