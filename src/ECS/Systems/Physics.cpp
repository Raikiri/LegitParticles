#include <string>
#include "Physics.h"
#include "../../Maths/AlgebraicMultigridSolver.h"
#include "../Context/MeshRendererData.h"
#include "../Context/InputData.h"
#include "../Context/CameraData.h"
#include "../Context/WindowData.h"
#include "SparsePhysics.h"
#include "../../Maths/StrainConstraints/StrainConstraints.h"
//#include "../../Maths/Tensor2.h"

namespace almost
{
  glm::vec4 Color32(glm::uint32_t col)
  {
    return glm::vec4(
      ((col & 0x000000ff) >> (0 * 8)) / 255.0f,
      ((col & 0x0000ff00) >> (1 * 8)) / 255.0f,
      ((col & 0x00ff0000) >> (2 * 8)) / 255.0f,
      ((col & 0xff000000) >> (3 * 8)) / 255.0f);
  }

  glm::vec4 LinearColor32(glm::uint32_t col)
  {
    return glm::pow(Color32(col), glm::vec4(2.2f));
  }

  PhysicsData InitPhysicsData()
  {
    PhysicsData physicsData;
    return physicsData;
  }

  void SubmitLine(almost::MeshRendererData& meshRendererData, glm::vec2 point0, glm::vec2 point1, float width, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    glm::vec2 dir = glm::normalize(point1 - point0);
    glm::vec2 tangent = glm::vec2(-dir.y, dir.x);

    meshRendererData.vertexData[vertexOffset + 0] = { glm::vec3(point0 + tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 1] = { glm::vec3(point0 - tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 2] = { glm::vec3(point1 - tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 3] = { glm::vec3(point1 + tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.indexData[indexOffset + 0] = glm::uint32_t(vertexOffset + 0);
    meshRendererData.indexData[indexOffset + 1] = glm::uint32_t(vertexOffset + 1);
    meshRendererData.indexData[indexOffset + 2] = glm::uint32_t(vertexOffset + 2);

    meshRendererData.indexData[indexOffset + 3] = glm::uint32_t(vertexOffset + 0);
    meshRendererData.indexData[indexOffset + 4] = glm::uint32_t(vertexOffset + 2);
    meshRendererData.indexData[indexOffset + 5] = glm::uint32_t(vertexOffset + 3);
    meshRendererData.verticesCount += 4;
    meshRendererData.indicesCount += 6;
  }

  void SubmitCircle(almost::MeshRendererData& meshRendererData, glm::vec2 center, float radius, size_t sectorsCount, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    for (size_t sectorIndex = 0; sectorIndex < sectorsCount; sectorIndex++)
    {
      float ang = (float(sectorIndex) / float(sectorsCount)) * 2.0f * glm::pi<float>();
      meshRendererData.vertexData[vertexOffset + sectorIndex] = { glm::vec3(center + glm::vec2(cos(ang), sin(ang)) * radius, 0.0f), glm::vec2(0.0f, 0.0f), color };

      if (sectorIndex < sectorsCount - 2)
        meshRendererData.indexData[indexOffset + sectorIndex * 3 + 0] = glm::uint32_t(vertexOffset);
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 1] = glm::uint32_t(vertexOffset + sectorIndex + 1);
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 2] = glm::uint32_t(vertexOffset + sectorIndex + 2);
    }
    meshRendererData.verticesCount += sectorsCount;
    meshRendererData.indicesCount += (sectorsCount - 2) * 3;

  }

  void SubmitTriangle(almost::MeshRendererData& meshRendererData, std::array<glm::vec2, 3> points, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    for (size_t vertexIndex = 0; vertexIndex < 3; vertexIndex++)
    {
      meshRendererData.vertexData[vertexOffset + vertexIndex] = { glm::vec3(points[vertexIndex], 0.0f), glm::vec2(0.0f, 0.0f), color };

      meshRendererData.indexData[indexOffset + vertexIndex] = glm::uint32_t(vertexOffset + vertexIndex);
    }
    meshRendererData.verticesCount += 3;
    meshRendererData.indicesCount += 3;
  }

  void ProcessPhysicsControls(
    WindowData& windowData,
    ParticleGroup particles,
    InputData& inputData,
    CameraData& cameraData)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();
    DefPosComponent* defPosComponents = particles.raw<DefPosComponent>();

    if (glfwGetKey(windowData.window, GLFW_KEY_SPACE))
    {
      for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
      {
        DefPosComponent& defPosComponent = defPosComponents[particleIndex];
        ParticleComponent& particleComponent = particleComponents[particleIndex];
        particleComponent.pos = defPosComponent.defPos;
#if defined(POSITION_BASED)
        particleComponent.prevPos = particleComponent.pos;
#else
        particleComponent.velocity = glm::vec2(0.0f);
#endif
      }
    }

    glm::vec2 camRightVec = glm::vec2(cos(cameraData.ang), sin(cameraData.ang));
    glm::vec2 camUpVec = glm::vec2(-camRightVec.y, camRightVec.x);
    bool isControlled = glfwGetMouseButton(windowData.window, GLFW_MOUSE_BUTTON_1);
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      DefPosComponent& defPosComponent = defPosComponents[particleIndex];
      MassComponent& massComponent = massComponents[particleIndex];
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      if (defPosComponent.isDraggable)
      {
        glm::vec2 dstPos = inputData.worldMousePos + camRightVec * defPosComponent.defPos.x + camUpVec * defPosComponent.defPos.y * 0.0f;
#if defined(POSITION_BASED)
        particleComponent.pos = isControlled ? dstPos : particleComponent.pos;
#else
        particleComponent.velocity = isControlled ? ((dstPos - particleComponent.pos) * 100.0f) : glm::vec2(0.0f);
#endif
      }
    }
  }

  void SolveLinksNonlinearGauss(ParticleGroup particles, LinkGroup links, legit::CpuProfiler& profiler)
  {
    auto physicsTask = profiler.StartScopedTask("[Physics] Nonlinear", legit::Colors::greenSea);

    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    LinkComponent* linkComponents = links.raw<LinkComponent>();

    for (int i = 0; i < 100; i++)
    {
      for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
      {
        size_t particleIndex0 = linkIndexComponents[linkIndex].indices[0];
        MassComponent& massComponent0 = massComponents[particleIndex0];
        ParticleComponent& particleComponent0 = particleComponents[particleIndex0];

        size_t particleIndex1 = linkIndexComponents[linkIndex].indices[1];
        MassComponent& massComponent1 = massComponents[particleIndex1];
        ParticleComponent& particleComponent1 = particleComponents[particleIndex1];

        glm::vec2 diff = particleComponent0.pos - particleComponent1.pos;
        glm::vec2 delta = diff / (glm::length(diff) + 1e-7f);
        float compInvMass = 1.0f / (massComponent0.invMass + massComponent1.invMass + 1e-7f);

        float deltaAcceleration = glm::dot(delta, particleComponent0.acceleration - particleComponent1.acceleration);
        //float deltaVelocity = glm::dot(delta, particleComponent0.velocity - particleComponent1.velocity);
        float deltaPos = glm::dot(delta, particleComponent0.pos - particleComponent1.pos) - linkComponents[linkIndex].defLength;

        float lambdaAcceleration = deltaAcceleration * compInvMass;
        //float lambdaVelocity = deltaVelocity * compInvMass;
        float lambdaPos = deltaPos * compInvMass;

        particleComponent0.pos += -delta * lambdaPos * massComponent0.invMass;
        particleComponent1.pos += delta * lambdaPos * massComponent1.invMass;
        /*particleComponent0.velocity += -delta * lambdaVelocity * massComponent0.invMass;
        particleComponent1.velocity += delta * lambdaVelocity * massComponent1.invMass;
        particleComponent0.acceleration += -delta * lambdaAcceleration * massComponent0.invMass;
        particleComponent1.acceleration += delta * lambdaAcceleration * massComponent1.invMass;*/
      }
    }
  }


  void SolveTriangles(ParticleGroup particles, TriangleGroup triangles, legit::CpuProfiler& profiler)
  {
    auto physicsTask = profiler.StartScopedTask("[Physics] Triangles", legit::Colors::greenSea);

    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    TriangleIndexComponent* triangleIndexComponents = triangles.raw<TriangleIndexComponent>();
    TriangleComponent* triangleComponents = triangles.raw<TriangleComponent>();

    for (int i = 0; i < 100; i++)
    {
      for (size_t triangleIndex = 0; triangleIndex < triangles.size(); triangleIndex++)
      {
        almost::StaticTensor<glm::vec2, almost::StrainDynamics::ParticleDim> worldPositions;
        almost::StaticTensor<float, almost::StrainDynamics::ParticleDim> invMasses;

        for (int p = 0; p < 1; p++)
        {
          for (auto m : almost::StrainDynamics::ParticleDim())
          {
            size_t particleIndex = triangleIndexComponents[triangleIndex].indices[m.value];
            worldPositions.Get(m) = particleComponents[particleIndex].pos;
            invMasses.Get(m) = massComponents[particleIndex].invMass;
          }
          almost::StrainDynamics::ProjectTriangle(worldPositions, invMasses, triangleComponents[triangleIndex].uvFromRef);
          for (auto m : almost::StrainDynamics::ParticleDim())
          {
            size_t particleIndex = triangleIndexComponents[triangleIndex].indices[m.value];
            particleComponents[particleIndex].pos = worldPositions.Get(m);
          }
        }
      }
    }
  }


  void ProcessPhysics(
    ParticleGroup particles,
    LinkGroup links,
    TriangleGroup triangles,
    almost::PhysicsData& physicsData,
    legit::CpuProfiler &profiler)
  {
    float dt = 1e-2f;
    glm::vec3 gravity = { 0.0f, -1000.0f, 0.0f };
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    ParticleIndexComponent* particleIndicesComponents = particles.raw<ParticleIndexComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    {
      auto physicsTask = profiler.StartScopedTask("[Physics] Links", legit::Colors::sunFlower);

      for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
      {
        MassComponent& massComponent = massComponents[particleIndex];
        ParticleComponent& particleComponent = particleComponents[particleIndex];
        particleComponent.acceleration = massComponent.usesGravity ? gravity : glm::vec3(0.0f, 0.0f, 0.0f);
        particleIndicesComponents[particleIndex].index = particleIndex;
      }

      LinkComponent* linkComponents = links.raw<LinkComponent>();
      LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
      for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
      {
        auto particleEntity0 = linkComponents[linkIndex].entities[0];
        auto particleEntity1 = linkComponents[linkIndex].entities[1];
        assert(particles.contains(particleEntity0));
        assert(particles.contains(particleEntity1));
        linkIndexComponents[linkIndex].indices[0] = particles.get<ParticleIndexComponent>(particleEntity0).index;
        linkIndexComponents[linkIndex].indices[1] = particles.get<ParticleIndexComponent>(particleEntity1).index;
      }

      TriangleComponent* triangleComponents = triangles.raw<TriangleComponent>();
      TriangleIndexComponent* triangleIndexComponents = triangles.raw<TriangleIndexComponent>();
      for (size_t triangleIndex = 0; triangleIndex < triangles.size(); triangleIndex++)
      {
        for (size_t particleNumber = 0; particleNumber < 3; particleNumber++)
        {
          auto particleEntity = triangleComponents[triangleIndex].entities[particleNumber];
          assert(particles.contains(particleEntity));
          triangleIndexComponents[triangleIndex].indices[particleNumber] = particles.get<ParticleIndexComponent>(particleEntity).index;
        }
      }
    }
    {

      /*for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
      {
        ParticleComponent& particleComponent = particleComponents[particleIndex];
        particleComponent.velocity += particleComponent.acceleration * dt;
      }*/
      //SolveLinksNonlinearGauss(particles, links, profiler);
      SolveTriangles(particles, triangles, profiler);

      //SolveLinksMultigrid(particles, links, profiler);
      //SolveLinksImplicitMultigrid(particles, links, dt, profiler);
      //SolveCollisions(particles, profiler);
    }

    {
      //auto physicsTask = profiler.StartScopedTask("[Physics] Integration", legit::Colors::greenSea);
      for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
      {
        ParticleComponent& particleComponent = particleComponents[particleIndex];
        #if defined(POSITION_BASED)
          glm::vec2 tmpPos = particleComponent.pos;
          particleComponent.pos += (particleComponent.pos - particleComponent.prevPos) + particleComponent.acceleration * dt * dt;
          particleComponent.prevPos = tmpPos;
        #else
          //particleComponent.velocity += particleComponent.acceleration * dt;
          particleComponent.pos += particleComponent.velocity * dt;
        #endif
      }
    }
  }

  void SubmitParticles(
    ParticleGroup &particles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    ParticleIndexComponent* particlesIndicesComponents = particles.raw<ParticleIndexComponent>();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      SubmitCircle(meshRendererData, particleComponents[particleIndex].pos, 3.0f, 15, LinearColor32(legit::Colors::belizeHole)/*favBlue*/);
    }
  }
  void SubmitLinks(
    ParticleGroup &particles,
    LinkGroup& links,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData)
  {
    LinkComponent* linkComponents = links.raw<LinkComponent>();
    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
    {
      auto particleEntity0 = linkComponents[linkIndex].entities[0];
      auto particleEntity1 = linkComponents[linkIndex].entities[1];
      assert(particles.contains(particleEntity0));
      assert(particles.contains(particleEntity1));
      ParticleComponent particleComponent0 = particles.get<ParticleComponent>(particleEntity0);
      ParticleComponent particleComponent1 = particles.get<ParticleComponent>(particleEntity1);
      SubmitLine(meshRendererData, particleComponent0.pos, particleComponent1.pos, 2.0f, LinearColor32(legit::Colors::orange)/*favOrange*/);
    }
  }

  void SubmitTriangles(
    ParticleGroup& particles,
    TriangleGroup& triangles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData)
  {
    TriangleComponent* triangleComponents = triangles.raw<TriangleComponent>();
    TriangleIndexComponent* triangleIndexComponents = triangles.raw<TriangleIndexComponent>();
    for (size_t triangleIndex = 0; triangleIndex < triangles.size(); triangleIndex++)
    {
      const auto& triangle = triangleComponents[triangleIndex];
      ParticleComponent particleComponent0 = particles.get<ParticleComponent>(triangle.entities[0]);
      ParticleComponent particleComponent1 = particles.get<ParticleComponent>(triangle.entities[1]);
      ParticleComponent particleComponent2 = particles.get<ParticleComponent>(triangle.entities[2]);

      glm::vec2 massCenter = glm::vec2(0);
      std::array<glm::vec2, 3> positions;
      for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
      {
        positions[vertexNumber] = particles.get<ParticleComponent>(triangle.entities[vertexNumber]).pos;
        massCenter += positions[vertexNumber] / 3.0f;
      }
      float ratio = 0.3f;
      for (size_t vertexNumber = 0; vertexNumber < 3; vertexNumber++)
      {
        positions[vertexNumber] = mix(positions[vertexNumber], massCenter, ratio);
      }
      SubmitTriangle(meshRendererData, positions, glm::vec4(glm::vec3(0.8f), 0.1f));
    }
  }
}