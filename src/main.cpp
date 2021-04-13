#include <entity\registry.hpp>

//#include <entt.hpp>


#include "Utils\ScopedCtx.h"
#include "ECS\Context\WindowData.h"
#include "ECS\Context\InputData.h"
#include "ECS\Context\RendererData.h"
#include "ECS\Context\MeshRendererData.h"
#include "ECS\Context\FullscreenRendererData.h"
#include "ECS\Context\PhysicsData.h"
#include "ECS\Context\CameraData.h"

#include "ECS\Systems\Window.h"
#include "ECS\Systems\Input.h"
#include "ECS\Systems\Renderer.h"
#include "ECS\Systems\MeshRenderer.h"
#include "ECS\Systems\FullscreenRenderer.h"
#include "ECS\Systems\Physics.h"
#include "ECS\Systems\Camera.h"

#include "ECS\Components\ParticleComponent.h"
#include "ECS\Components\LinkComponent.h"
#include "ECS\Components\TriangleComponent.h"

#include "ECS\Scenes\ParticlesScene.h"

int main(int argsCount, char **args)
{
  entt::registry reg;
  auto windowHandle = ScopedCtx<almost::WindowData>(reg, almost::InitWindow());
  auto rendererHandle = ScopedCtx<almost::RendererData>(reg, almost::InitRenderer(reg.ctx<almost::WindowData>()));
  auto inputHandle = ScopedCtx<almost::InputData>(reg, almost::InitInputData(reg.ctx<almost::WindowData>()));
  auto cameraHandle = ScopedCtx<almost::CameraData>(reg, almost::InitCameraData());
  auto meshRendererHandle = ScopedCtx<almost::MeshRendererData>(reg, almost::InitMeshRendererData(reg.ctx<almost::RendererData>()));
  auto physicsHandle = ScopedCtx<almost::PhysicsData>(reg, almost::InitPhysicsData());
  auto fullscreenRendererHandle = ScopedCtx<almost::FullscreenRendererData>(reg, almost::InitFullscreenRendererData(reg.ctx<almost::RendererData>()));
  almost::MapMeshRenderer(reg.ctx<almost::MeshRendererData>());

  CreateGround(reg);
  CreateClothPhysicsMesh(reg);

  while (!reg.ctx<almost::InputData>().isWindowClosed)
  {
    auto& windowData = reg.ctx<almost::WindowData>();
    auto& rendererData = reg.ctx<almost::RendererData>();
    auto& fullscreenRendererData = reg.ctx<almost::FullscreenRendererData>();
    auto& meshRendererData = reg.ctx<almost::MeshRendererData>();
    auto& cameraData = reg.ctx<almost::CameraData>();
    auto& physicsData = reg.ctx<almost::PhysicsData>();
    auto& inputData = reg.ctx<almost::InputData>();

    //almost::ProcessPhysicsData(reg.ctx<almost::PhysicsData>(), reg.ctx<almost::MeshRendererData>());

    auto particleGroup = reg.group<almost::ParticleComponent, almost::ParticleIndexComponent, almost::MassComponent, almost::DefPosComponent>();
    auto linkGroup = reg.group<almost::LinkComponent, almost::LinkIndexComponent>();
    auto triangleGroup = reg.group<almost::TriangleComponent, almost::TriangleIndexComponent>();

    almost::StartFrame(windowData, rendererData);
    reg.ctx<almost::MeshRendererData>().verticesCount = 0;
    reg.ctx<almost::MeshRendererData>().indicesCount = 0;

    {
      auto inputTask = rendererData.inFlightQueue->GetCpuProfiler().StartScopedTask("Input", legit::Colors::greenSea);
      almost::ProcessInput(reg.ctx<almost::WindowData>(), reg.ctx<almost::InputData>());
      almost::UpdateCamera(reg.ctx<almost::WindowData>(), reg.ctx<almost::RendererData>(), reg.ctx<almost::InputData>(), reg.ctx<almost::CameraData>());
      almost::ProcessPhysicsControls(windowData, particleGroup, inputData, cameraData);
    }

    {
      almost::ProcessPhysics(particleGroup, linkGroup, triangleGroup, physicsData, rendererData.inFlightQueue->GetCpuProfiler());
    }

    {
      auto submitTask = rendererData.inFlightQueue->GetCpuProfiler().StartScopedTask("Submit", legit::Colors::turqoise);

      almost::SubmitTriangles(particleGroup, triangleGroup, physicsData, meshRendererData);
      almost::SubmitLinks(particleGroup, linkGroup, physicsData, meshRendererData);
      almost::SubmitParticles(particleGroup, physicsData, meshRendererData);
      almost::RenderFullscreen(windowData, fullscreenRendererData, rendererData, cameraData);
      almost::RenderMeshes(windowData, meshRendererData, rendererData, cameraData);
    }
    almost::FinishFrame(rendererData);
  }
  almost::UnmapMeshRenderer(reg.ctx<almost::MeshRendererData>());
  reg.ctx<almost::RendererData>().core->WaitIdle();
  glfwDestroyWindow(reg.ctx<almost::WindowData>().window);
  glfwTerminate();

  return 0;
}