#include <string>
#include "PhysicsCommon.h"
#include "../Context/MeshRendererData.h"
#include "../Context/InputData.h"
#include "../Context/CameraData.h"
#include "../Context/WindowData.h"
#include "imgui.h"

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



  void ProcessPhysicsControls(
    WindowData& windowData,
    almost::ParticleGroup::Type particles,
    InputData& inputData,
    CameraData& cameraData)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();
    DefPosComponent* defPosComponents = particles.raw<DefPosComponent>();

    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_SPACE))
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
    bool isControlled = ImGui::IsMouseDown(0) && !ImGui::IsAnyWindowFocused();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      DefPosComponent& defPosComponent = defPosComponents[particleIndex];
      MassComponent& massComponent = massComponents[particleIndex];
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      if (defPosComponent.isDraggable)
      {
        glm::vec2 dstPos = inputData.worldMousePos + camRightVec * defPosComponent.defPos.x + camUpVec * defPosComponent.defPos.y * 0.0f;
        particleComponent.pos = isControlled ? dstPos : particleComponent.pos;
#if defined(POSITION_BASED)
        //particleComponent.pos = isControlled ? dstPos : particleComponent.pos;
#else
        //particleComponent.velocity = isControlled ? ((dstPos - particleComponent.pos) * 100.0f) : glm::vec2(0.0f);
#endif
      }
    }
  }

  void PreStep(
    ParticleGroup::Type particleGroup,
    LinkGroup::Type linkGroup,
    TriangleGroup::Type triangleGroup)
  {
    glm::vec3 gravity = { 0.0f, -100.0, 0.0f };
    ParticleComponent* particleComponents = particleGroup.raw<ParticleComponent>();
    ParticleIndexComponent* particleIndicesComponents = particleGroup.raw<ParticleIndexComponent>();
    MassComponent* massComponents = particleGroup.raw<MassComponent>();

    {
      for (size_t particleIndex = 0; particleIndex < particleGroup.size(); particleIndex++)
      {
        MassComponent& massComponent = massComponents[particleIndex];
        ParticleComponent& particleComponent = particleComponents[particleIndex];
        particleComponent.acceleration = massComponent.usesGravity ? gravity : glm::vec3(0.0f, 0.0f, 0.0f);
        particleIndicesComponents[particleIndex].index = particleIndex;
      }

      LinkComponent* linkComponents = linkGroup.raw<LinkComponent>();
      LinkIndexComponent* linkIndexComponents = linkGroup.raw<LinkIndexComponent>();
      for (size_t linkIndex = 0; linkIndex < linkGroup.size(); linkIndex++)
      {
        auto particleEntity0 = linkComponents[linkIndex].entities[0];
        auto particleEntity1 = linkComponents[linkIndex].entities[1];
        assert(particleGroup.contains(particleEntity0));
        assert(particleGroup.contains(particleEntity1));
        linkIndexComponents[linkIndex].indices[0] = particleGroup.get<ParticleIndexComponent>(particleEntity0).index;
        linkIndexComponents[linkIndex].indices[1] = particleGroup.get<ParticleIndexComponent>(particleEntity1).index;
      }

      TriangleComponent* triangleComponents = triangleGroup.raw<TriangleComponent>();
      TriangleIndexComponent* triangleIndexComponents = triangleGroup.raw<TriangleIndexComponent>();
      for (size_t triangleIndex = 0; triangleIndex < triangleGroup.size(); triangleIndex++)
      {
        for (size_t particleNumber = 0; particleNumber < 3; particleNumber++)
        {
          auto particleEntity = triangleComponents[triangleIndex].entities[particleNumber];
          assert(particleGroup.contains(particleEntity));
          triangleIndexComponents[triangleIndex].indices[particleNumber] = particleGroup.get<ParticleIndexComponent>(particleEntity).index;
        }
      }
    }
  }
  
  void ProjectLinks(
    ParticleComponent* particleComponents, MassComponent* massComponents, size_t particlesCount,
    LinkComponent* linkComponents, LinkIndexComponent* linkIndexComponents, size_t linksCount)
  {
    for (size_t linkIndex = 0; linkIndex < linksCount; linkIndex++)
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
      //if (deltaPos < 0) continue;
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
  
  void IntegrateParticles(ParticleComponent* particleComponents, size_t particlesCount, float dt)
  {
    //glm::vec2 acc = {0.0f, 10.0f};
    for (size_t particleIndex = 0; particleIndex < particlesCount; particleIndex++)
    {
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      #if defined(POSITION_BASED)
        glm::vec2 tmpPos = particleComponent.pos;
        particleComponent.pos += (particleComponent.pos - particleComponent.prevPos) + particleComponent.acceleration * dt * dt;
        particleComponent.prevPos = tmpPos;
      #else
        particleComponent.velocity += particleComponent.acceleration * dt;
        particleComponent.pos += particleComponent.velocity * dt;
      #endif
    }
  }
  
  void SubmitParticles(
    ParticleGroup::Type particles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      SubmitCircle(meshRendererData, particleComponents[particleIndex].pos, 3.0f, 15, LinearColor32(legit::Colors::belizeHole)/*favBlue*/);
    }
  }
  void SubmitLinks(
    ParticleGroup::Type particles,
    LinkGroup::Type links,
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
    ParticleGroup::Type particles,
    TriangleGroup::Type triangles,
    almost::PhysicsData &physicsData,
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