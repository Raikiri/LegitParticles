#pragma once
#include "../../LegitProfiler/ProfilerTask.h"
#include "../../LegitVulkan/CpuProfiler.h"

#include "../Context/PhysicsData.h"
//#include "../Context/MeshRendererData.h"
#include "../Components/ParticleComponent.h"
#include "../Components/LinkComponent.h"
#include "../Components/TriangleComponent.h"
#include <entity/registry.hpp>
#include "../../Utils/GroupArg.h"
#include "../../Maths/StackStorage.h"
namespace almost
{
  const glm::vec4 favBlue = glm::vec4(glm::pow(glm::vec3(0.0, 0.5, 0.75), glm::vec3(2.2f)), 1.0f);
  const glm::vec4 favOrange = glm::vec4(glm::pow(glm::vec3(0.8f, 0.4f, 0.1f), glm::vec3(2.2f)), 1.0f);


  PhysicsData InitPhysicsData();
  struct WindowData;
  struct InputData;
  struct CameraData;
  struct MeshRendererData;

  template<typename ...Components>
  struct Group
  {
    using Type = entt::group_type<Components...>;
    static Type Get(entt::registry &reg)
    {
      return reg.group<Components...>();
    }
  };
  using ParticleGroup = Group<almost::ParticleComponent, almost::ParticleIndexComponent, almost::MassComponent, almost::DefPosComponent, almost::CoarseMultigridComponent, almost::FineMultigridComponent>;
  using LinkGroup = Group<almost::LinkComponent, almost::LinkIndexComponent>;
  using TriangleGroup = Group<almost::TriangleComponent, almost::TriangleIndexComponent>;

  void ProcessPhysicsControls(
    WindowData& windowData,
    ParticleGroup::Type particles,
    InputData& inputData,
    CameraData& cameraData);
  

  void ProcessPhysics(
    std::vector<entt::registry> &regLayers,
    legit::CpuProfiler& profiler);

  void SubmitParticles(
    ParticleGroup::Type particles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData);
  void SubmitLinks(
    ParticleGroup::Type particles,
    LinkGroup::Type links,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData);
  void SubmitTriangles(
    ParticleGroup::Type particles,
    TriangleGroup::Type triangles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData);
}