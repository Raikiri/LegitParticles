#pragma once
#include "../../LegitVulkan/ProfilerTask.h"
#include "../../LegitVulkan/CpuProfiler.h"

#include "../Context/PhysicsData.h"
//#include "../Context/MeshRendererData.h"
#include "../Components/ParticleComponent.h"
#include "../Components/LinkComponent.h"
#include "../Components/TriangleComponent.h"
#include <entity\registry.hpp>
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

  using ParticleGroup = entt::group_type<almost::ParticleComponent, almost::ParticleIndexComponent, almost::MassComponent, almost::DefPosComponent>;
  using LinkGroup = entt::group_type<almost::LinkComponent, almost::LinkIndexComponent>;
  using TriangleGroup = entt::group_type<almost::TriangleComponent, almost::TriangleIndexComponent>;

  void ProcessPhysicsControls(
    WindowData& windowData,
    ParticleGroup particles,
    InputData& inputData,
    CameraData& cameraData);
  

  void ProcessPhysics(
    ParticleGroup particles,
    LinkGroup links,
    TriangleGroup triangles,
    almost::PhysicsData& physicsData,
    legit::CpuProfiler& profiler);

  void SubmitParticles(
    ParticleGroup& particles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData);
  void SubmitLinks(
    ParticleGroup& particles,
    LinkGroup& links,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData);
  void SubmitTriangles(
    ParticleGroup& particles,
    TriangleGroup& triangles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData);
}