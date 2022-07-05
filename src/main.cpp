#include <entity\registry.hpp>
#include "LegitVulkan/LegitVulkan.h"
#include "LegitVulkan/CoreImpl.h"

//#include <entt.hpp>


#include "Utils\ScopedCtx.h"
#include "ECS\Context\WindowData.h"
#include "ECS\Context\InputData.h"
#include "ECS\Context\RendererData.h"
#include "ECS\Context\MeshRendererData.h"
#include "ECS\Context\FullscreenRendererData.h"
#include "ECS\Context\ShaderTestbedRendererData.h"
#include "ECS\Context\PhysicsData.h"
#include "ECS\Context\AnimationData.h"
#include "ECS\Context\CameraData.h"
#include "ECS\Context\PhysicsAnimationData.h"
#include "ECS\Systems\Window.h"
#include "ECS\Systems\Input.h"
#include "ECS\Systems\Renderer.h"
#include "ECS\Systems\MeshRenderer.h"
#include "ECS\Systems\FullscreenRenderer.h"
#include "ECS\Systems\ShaderTestbedRenderer.h"
#include "ECS\Systems\Physics.h"
#include "ECS\Systems\SkeletalAnimation.h"
#include "ECS\Systems\PhysicsAnimation.h"
#include "ECS\Systems\Camera.h"

#include "ECS\Components\ParticleComponent.h"
#include "ECS\Components\LinkComponent.h"
#include "ECS\Components\TriangleComponent.h"

#include "ECS\Scenes\MultigridClothScene.h"


int main(int argsCount, char **args)
{
  entt::registry reg;
  auto windowHandle = ScopedCtx<almost::WindowData>(reg, almost::InitWindow());
  auto rendererHandle = ScopedCtx<almost::RendererData>(reg, almost::InitRenderer(reg.ctx<almost::WindowData>()));
  auto inputHandle = ScopedCtx<almost::InputData>(reg, almost::InitInputData(reg.ctx<almost::WindowData>()));
  auto cameraHandle = ScopedCtx<almost::CameraData>(reg, almost::InitCameraData());
  auto meshRendererHandle = ScopedCtx<almost::MeshRendererData>(reg, almost::InitMeshRendererData(reg.ctx<almost::RendererData>()));
  auto physicsHandle = ScopedCtx<almost::PhysicsData>(reg, almost::InitPhysicsData());
  auto animationHandle = ScopedCtx<almost::AnimationData>(reg, almost::InitAnimationData());
  auto physicsAnimationHandle = ScopedCtx<almost::PhysicsAnimationData>(reg, almost::InitPhysicsAnimationData());
  auto fullscreenRendererHandle = ScopedCtx<almost::FullscreenRendererData>(reg, almost::InitFullscreenRendererData(reg.ctx<almost::RendererData>()));
  auto tesetbedRendererHandle = ScopedCtx<almost::ShaderTestbedRendererData>(reg, almost::InitShaderTestbedRendererData(reg.ctx<almost::RendererData>()));
  almost::MapMeshRenderer(reg.ctx<almost::MeshRendererData>());

  std::vector<entt::registry> regLayers;
  regLayers.resize(6);

  CreateGround(regLayers[0]);
  //CreateClothPhysicsMesh(reg);
  //CreateMultigridPhysicsMesh(regLayers, { -200, -200 }, { 200, 200 }, { 65, 65 });

  while (!reg.ctx<almost::InputData>().isWindowClosed)
  {
    auto& windowData = reg.ctx<almost::WindowData>();
    auto& rendererData = reg.ctx<almost::RendererData>();
    auto& fullscreenRendererData = reg.ctx<almost::FullscreenRendererData>();
    auto& shaderTestbedRendererData = reg.ctx<almost::ShaderTestbedRendererData>();
    auto& meshRendererData = reg.ctx<almost::MeshRendererData>();
    auto& cameraData = reg.ctx<almost::CameraData>();
    auto& physicsData = reg.ctx<almost::PhysicsData>();
    auto& animationData = reg.ctx<almost::AnimationData>();
    auto& physicsAnimationData = reg.ctx<almost::PhysicsAnimationData>();
    auto& inputData = reg.ctx<almost::InputData>();

    //almost::ProcessPhysicsData(reg.ctx<almost::PhysicsData>(), reg.ctx<almost::MeshRendererData>());


    almost::StartFrame(windowData, rendererData);
    reg.ctx<almost::MeshRendererData>().verticesCount = 0;
    reg.ctx<almost::MeshRendererData>().indicesCount = 0;

    {
      auto inputTask = rendererData.inFlightQueue->GetCpuProfiler().StartScopedTask("Input", legit::Colors::greenSea);
      almost::ProcessInput(reg.ctx<almost::WindowData>(), reg.ctx<almost::InputData>());
      almost::UpdateCamera(reg.ctx<almost::WindowData>(), reg.ctx<almost::RendererData>(), reg.ctx<almost::InputData>(), reg.ctx<almost::CameraData>());
      for (auto& layer : regLayers)
      {
        almost::ProcessPhysicsControls(windowData, almost::ParticleGroup::Get(layer), inputData, cameraData);
      }
    }

    {
      almost::ProcessPhysics(regLayers, rendererData.inFlightQueue->GetCpuProfiler());
    }
    {
      almost::ProcessPhysicsAnimation(physicsAnimationData, inputData, windowData);
    }

    {
      auto submitTask = rendererData.inFlightQueue->GetCpuProfiler().StartScopedTask("Submit", legit::Colors::turqoise);

      for (auto& layer : regLayers)
      {
        almost::SubmitTriangles(almost::ParticleGroup::Get(layer), almost::TriangleGroup::Get(layer), physicsData, meshRendererData);
        almost::SubmitLinks(almost::ParticleGroup::Get(layer), almost::LinkGroup::Get(layer), physicsData, meshRendererData);
        almost::SubmitParticles(almost::ParticleGroup::Get(layer), physicsData, meshRendererData);
      }

      //almost::SubmitAnimation(animationData, meshRendererData);
      almost::SubmitPhysicsAnimation(physicsAnimationData, inputData, meshRendererData);
      almost::RenderFullscreen(windowData, fullscreenRendererData, rendererData, cameraData);
      almost::RenderMeshes(windowData, meshRendererData, rendererData, cameraData);

      //almost::RenderShaderTestbed(windowData, shaderTestbedRendererData, rendererData, cameraData);
    }
    almost::FinishFrame(rendererData);
  }
  almost::UnmapMeshRenderer(reg.ctx<almost::MeshRendererData>());
  reg.ctx<almost::RendererData>().core->WaitIdle();
  glfwDestroyWindow(reg.ctx<almost::WindowData>().window);
  glfwTerminate();

  return 0;
}