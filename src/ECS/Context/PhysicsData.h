#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../Maths/SparseMatrix.h"

namespace almost
{
  struct PhysicsData
  {
    using LinkIndex = size_t;
    using ParticleIndex = size_t;

    using JacobianMatrix = almost::SparseMatrix<LinkIndex, ParticleIndex, glm::vec2>;
    using JacobianTransposedMatrix = almost::SparseMatrix<ParticleIndex, LinkIndex, glm::vec2>;
    using MassMatrix = almost::SparseMatrix<ParticleIndex, ParticleIndex, float>;
    using DeltaMatrix = almost::SparseMatrix<ParticleIndex, char, glm::vec2>;

    using SystemMatrix = almost::SparseMatrix<LinkIndex, LinkIndex, float>;
    using RightSideMatrix = almost::SparseMatrix<LinkIndex, char, float>;
  };
}