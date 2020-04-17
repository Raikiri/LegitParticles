#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../Maths/SparseMatrix.h"

namespace almost
{
  struct Vector2f
  {

    Vector2f& operator += (Vector2f& other)
    {
      vec += other.vec;
      return *this;
    }
    glm::vec2 vec;
  };
  struct PhysicsData
  {
    using LinkIndex = size_t;
    using ParticleIndex = size_t;

    using JacobianMatrix = almost::SparseMatrix<LinkIndex, ParticleIndex, Vector2f>;
    using JacobianTransposedMatrix = almost::SparseMatrix<ParticleIndex, LinkIndex, Vector2f>;
    using MassMatrix = almost::SparseMatrix<ParticleIndex, ParticleIndex, float>;
    using DeltaMatrix = almost::SparseMatrix<ParticleIndex, char, Vector2f>;

    using SystemMatrix = almost::SparseMatrix<LinkIndex, LinkIndex, float>;
    using RightSideMatrix = almost::SparseMatrix<LinkIndex, char, float>;
  };
}