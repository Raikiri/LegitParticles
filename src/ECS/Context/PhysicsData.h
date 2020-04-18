#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../Maths/SparseMatrix.h"

namespace almost
{
  struct PhysicsData
  {
    struct ParticleDimension
    {
      using IndexType = size_t;
      size_t size;
    };
    struct LinkDimension
    {
      using IndexType = size_t;
      size_t size;
    };
    struct OneDimension
    {
      using IndexType = char;
      const static size_t size = 1;
    };

    using JacobianMatrix = almost::SparseMatrix<LinkDimension, ParticleDimension, glm::vec2>;
    using JacobianTransposedMatrix = almost::SparseMatrix<ParticleDimension, LinkDimension, glm::vec2>;
    using MassMatrix = almost::SparseMatrix<ParticleDimension, ParticleDimension, float>;
    using DeltaMatrix = almost::SparseMatrix<ParticleDimension, OneDimension, glm::vec2>;

    using SystemMatrix = almost::SparseMatrix<LinkDimension, LinkDimension, float>;
    using RightSideMatrix = almost::SparseMatrix<LinkDimension, OneDimension, float>;
  };
}