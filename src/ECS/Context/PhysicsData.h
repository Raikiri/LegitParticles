#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../Maths/SparseMatrix.h"
#include "../../Maths/Tensor2.h"
namespace almost
{
  struct PhysicsData
  {
    struct ParticleDimension
    {
      using IndexType = size_t;
      size_t size;
    };
    struct JointDimension
    {
      using IndexType = size_t;
      size_t size;
    };
    struct OneDimension
    {
      using IndexType = char;
      const static size_t size = 1;
    };

    using JacobianMatrix = almost::SparseMatrix<JointDimension, ParticleDimension, glm::vec2>;
    using JacobianTransposedMatrix = almost::SparseMatrix<ParticleDimension, JointDimension, glm::vec2>;
    using MassMatrix = almost::SparseMatrix<ParticleDimension, ParticleDimension, float>;
    using DeltaMatrix = almost::SparseMatrix<ParticleDimension, OneDimension, glm::vec2>;

    using SystemMatrix = almost::SparseMatrix<JointDimension, JointDimension, float>;
    using RightSideMatrix = almost::SparseMatrix<JointDimension, OneDimension, float>;

    using ParticleTensorMatrix = almost::SparseMatrix<ParticleDimension, ParticleDimension, almost::Tensor2f>;
    using ParticleScalarMatrix = almost::SparseMatrix<ParticleDimension, ParticleDimension, float>;
  };
}