#pragma once
#include "../TensorMaths.h"
#include "../../Utils/GlmInclude.h"

namespace almost
{

  namespace StrainDynamics
  {
    struct WorldTag {};
    using WorldDim = almost::StaticDimension<WorldTag, 2>;

    struct RefTag {};
    using RefDim = almost::StaticDimension<RefTag, 2>;

    struct UvTag {};
    using UvDim = almost::StaticDimension<UvTag, 2>;

    struct ParticleTag {};
    using ParticleDim = almost::StaticDimension<ParticleTag, 3>;

    struct Basis {};
    using BasisDim = almost::DynamicDimension<Basis>;

    using GlmIndex = glm::length_t;

    almost::StaticTensor<float, UvDim, RefDim> BuildUvFromRefMatrix(almost::StaticTensor<glm::vec2, ParticleDim> refParticles);


    float GetStrain(
      const almost::StaticTensor<float, WorldDim, RefDim>& worldFromRef, //F
      almost::Index<RefDim> i, almost::Index<RefDim> j);

    const almost::StaticTensor<glm::vec2, ParticleDim> GetStrainParticleGradients(
      const almost::StaticTensor<float, WorldDim, RefDim>& worldFromRef, //F
      const almost::StaticTensor<float, UvDim, RefDim>& uvFromRef, //Q^-1
      almost::Index<RefDim> i, almost::Index<RefDim> j);

    struct TriangleJoint
    {
      almost::StaticTensor<glm::vec2, ParticleDim> jacobian;
      float rhs;
    };

    const TriangleJoint BuildSqrtStretchTriangleJoint(const almost::StaticTensor<glm::vec2, ParticleDim> strainGradients, float strain);

    const almost::StaticTensor<glm::vec2, ParticleDim> SolveTriangleJoint(almost::StaticTensor<float, ParticleDim> invMasses, const TriangleJoint& joint);

    struct TensorIndices
    {
      TensorIndices(size_t i, size_t j) :
        i(i),
        j(j)
      {

      }
      almost::Index<RefDim> i, j;
    };

    void ProjectTriangle(
      almost::StaticTensor<glm::vec2, ParticleDim>& worldParticles,
      const almost::StaticTensor<float, ParticleDim>& invMasses,
      const almost::StaticTensor<float, UvDim, RefDim> uvFromRef /*Q^-1*/);
  }
}