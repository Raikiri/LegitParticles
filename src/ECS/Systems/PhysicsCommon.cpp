#include <string>
#include "PhysicsCommon.h"
#include "../Context/MeshRendererData.h"
#include "../Context/InputData.h"
#include "../Context/CameraData.h"
#include "../Context/WindowData.h"

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
    bool isControlled = glfwGetMouseButton(windowData.window->glfw_window, GLFW_MOUSE_BUTTON_1);
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

  void SubmitParticles(
    ParticleGroup::Type particles,
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