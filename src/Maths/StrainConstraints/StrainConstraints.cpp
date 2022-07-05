#include "StrainConstraints.h"

namespace almost
{

  namespace StrainDynamics
  {
    almost::StaticTensor<float, UvDim, RefDim> BuildUvFromRefMatrix(almost::StaticTensor<glm::vec2, ParticleDim> refParticles)
    {
      glm::mat2 refFromUv;

      auto pos0 = refParticles.Get(almost::Index<ParticleDim>(0));
      for (auto n : RefDim())
        for (auto m : UvDim())
          refFromUv[GlmIndex(m.value)][GlmIndex(n.value)] = (refParticles.Get(almost::Index<ParticleDim>(m.value + 1)) - pos0)[GlmIndex(n.value)]; //Q

      auto uvFromRef = glm::inverse(refFromUv);

      auto uvFromRefTensor = almost::MakeStaticTensor<float, UvDim, RefDim>();

      for (auto n : RefDim())
        for (auto m : UvDim())
        uvFromRefTensor.Get(m, n) = uvFromRef[GlmIndex(n.value)][GlmIndex(m.value)];

      return uvFromRefTensor;
    }


    float GetStrain(
      const almost::StaticTensor<float, WorldDim, RefDim>& worldFromRef, //F
      almost::Index<RefDim> i, almost::Index<RefDim> j)
    {
      float strain = 0.0f;
      for (auto k : WorldDim())
      {
        strain += worldFromRef.Get(k, i) * worldFromRef.Get(k, j); //S = G + I
      }
      return strain;
    }

    const almost::StaticTensor<glm::vec2, ParticleDim> GetStrainParticleGradients(
      const almost::StaticTensor<float, WorldDim, RefDim>& worldFromRef, //F
      const almost::StaticTensor<float, UvDim, RefDim>& uvFromRef, //Q^-1
      almost::Index<RefDim> i, almost::Index<RefDim> j)
    {
      almost::StaticTensor<glm::vec2, ParticleDim> strainGradients;
      glm::vec2 sum{ 0 };

      glm::vec2 zeroGrad{ 0 };
      for (auto m : UvDim())
      {
        auto p = almost::Index<ParticleDim>(m.value + 1);

        for (auto n : WorldDim())
        {
          float val =
            uvFromRef.Get(m, i) * worldFromRef.Get(n, j) +
            uvFromRef.Get(m, j) * worldFromRef.Get(n, i);

          strainGradients.Get(p)[GlmIndex(n.value)] = val;
        }
        zeroGrad -= strainGradients.Get(p);
      }
      strainGradients.Get(almost::Index<ParticleDim>(0)) = zeroGrad;
      return strainGradients;
    }

    const TriangleJoint BuildSqrtStretchTriangleJoint(const almost::StaticTensor<glm::vec2, ParticleDim> strainGradients, float strain)
    {
      TriangleJoint joint;
      float strainSqrt = sqrt(strain);
      float mult = 1.0f / (2.0f * strainSqrt);
      for (auto m : ParticleDim())
      {
        joint.jacobian.Get(m) = strainGradients.Get(m) * mult;
      }
      joint.rhs = std::max(0.0f, strainSqrt - 1.0f);
      return joint;
    }

    const almost::StaticTensor<glm::vec2, ParticleDim> SolveTriangleJoint(almost::StaticTensor<float, ParticleDim> invMasses, const TriangleJoint& joint)
    {
      float compMass = 1e-7f;
      for (auto m : ParticleDim())
      {
        compMass += glm::dot(joint.jacobian.Get(m), joint.jacobian.Get(m)) * invMasses.Get(m);
      }
      float lambda = joint.rhs / compMass;

      almost::StaticTensor<glm::vec2, ParticleDim> particleOffsets;
      for (auto m : ParticleDim())
      {
        particleOffsets.Get(m) = -invMasses.Get(m) * lambda * joint.jacobian.Get(m);
      }
      return particleOffsets;
    }

    void ProjectTriangle(almost::StaticTensor<glm::vec2, ParticleDim>& worldParticles, const almost::StaticTensor<float, ParticleDim>& invMasses, const almost::StaticTensor<float, UvDim, RefDim> uvFromRef /*Q^-1*/)
    {


      auto worldFromUv = almost::MakeStaticTensor<float, WorldDim, UvDim>();

      auto pos0 = worldParticles.Get(almost::Index<ParticleDim>(0));
      for (auto n : WorldDim())
        for (auto m : UvDim())
          worldFromUv.Get(n, m) = worldParticles.Get(almost::Index<ParticleDim>(m.value + 1))[GlmIndex(n.value)] - pos0[GlmIndex(n.value)]; //P

      auto worldFromRef = almost::MakeStaticTensor<float, WorldDim, RefDim>();
      for (auto k : WorldDim())
      {
        for (auto i : RefDim())
        {
          worldFromRef.Get(k, i) = 0.0f;
          for (auto x : UvDim())
          {
            worldFromRef.Get(k, i) += worldFromUv.Get(k, x) * uvFromRef.Get(x, i); //F
          }
        }
      }

      auto stretchIndices = { TensorIndices{0, 0}, TensorIndices{1, 1} };
      auto shearIndices = { TensorIndices{0, 1} };
      for (auto index : stretchIndices)
      {
        float strain = GetStrain(worldFromRef, index.i, index.j);
        auto strainGradients = GetStrainParticleGradients(worldFromRef, uvFromRef, index.i, index.j);
        auto joint = BuildSqrtStretchTriangleJoint(strainGradients, strain);
        auto particleOffsets = SolveTriangleJoint(invMasses, joint);
        for (auto m : ParticleDim())
          worldParticles.Get(m) += particleOffsets.Get(m);
      }
    }
  }
}