#include "../Context/RendererData.h"
#include "../Context/FullscreenRendererData.h"
#include "../Context/CameraData.h"

#include <memory>
namespace almost
{

  void ReloadShader(almost::RendererData& rendererData, FullscreenRendererData& fullscreenRendererData)
  {
    auto* core = rendererData.core.get();
    vk::PhysicalDevice physicalDevice = core->GetPhysicalDevice();
    vk::Device logicalDevice = core->GetLogicalDevice();

    fullscreenRendererData.shader.vertex.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/screenspaceQuad.vert.spv"));
    fullscreenRendererData.shader.fragment.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/fullscreenRenderer.frag.spv"));
    fullscreenRendererData.shader.program.reset(new legit::ShaderProgram(fullscreenRendererData.shader.vertex.get(), fullscreenRendererData.shader.fragment.get()));
  }
  FullscreenRendererData InitFullscreenRendererData(almost::RendererData& rendererData)
  {

    FullscreenRendererData fullscreenRendererData;
    ReloadShader(rendererData, fullscreenRendererData);

    return fullscreenRendererData;
  }

  void RenderFullscreen(const almost::WindowData& windowData, almost::FullscreenRendererData& fullscreenRendererData, almost::RendererData& rendererData, const almost::CameraData &cameraData)
  {
    #pragma pack(push, 1)
    struct GlobalDataBuffer
    {
      glm::vec2 camPos;
      glm::vec2 camFov;
      glm::vec2 viewportSize;
      float camAngle;
    };
    #pragma pack(pop)
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_V))
    {
      ReloadShader(rendererData, fullscreenRendererData);
    }
    rendererData.core->GetRenderGraph()->AddPass(legit::RenderGraph::RenderPassDesc()
      .SetColorAttachments({ { rendererData.frameInfo.swapchainImageViewProxyId, vk::AttachmentLoadOp::eClear, rendererData.clearColor}  })
      .SetProfilerInfo(legit::Colors::turqoise, "Grid")
      .SetRenderAreaExtent(rendererData.inFlightQueue->GetImageSize())
      .SetRecordFunc( [&fullscreenRendererData, &rendererData, &cameraData](legit::RenderGraph::RenderPassContext passContext)
    {
      const static uint32_t ShaderDataSetIndex = 0;
      const static uint32_t DrawCallDataSetIndex = 1;

      auto shader = fullscreenRendererData.shader.program.get();

      auto pipeineInfo = rendererData.core->GetPipelineCache()->BindGraphicsPipeline(passContext.GetCommandBuffer(), passContext.GetRenderPass()->GetHandle(), legit::DepthSettings::Disabled(), { legit::BlendSettings::Opaque() }, legit::VertexDeclaration(), vk::PrimitiveTopology::eTriangleFan, shader);
      {
        const legit::DescriptorSetLayoutKey *shaderDataSetInfo = shader->GetSetInfo(ShaderDataSetIndex);
        auto shaderData = rendererData.frameInfo.memoryPool->BeginSet(shaderDataSetInfo);
        {
          auto shaderDataBuffer = rendererData.frameInfo.memoryPool->GetUniformBufferData<GlobalDataBuffer>("GlobalDataBuffer");

          shaderDataBuffer->camPos = cameraData.pos;
          shaderDataBuffer->camFov = cameraData.fov;
          auto viewportExtent = rendererData.inFlightQueue->GetImageSize();
          shaderDataBuffer->viewportSize = glm::vec2(viewportExtent.width, viewportExtent.height);
          shaderDataBuffer->camAngle = cameraData.ang;
        }
        rendererData.frameInfo.memoryPool->EndSet();

        /*std::vector<legit::ImageSamplerBinding> imageSamplerBindings;
        auto directLightImageView = passContext.GetImageView(this->viewportResources->directLight.imageViewProxy->Id());
        imageSamplerBindings.push_back(shaderDataSetInfo->MakeImageSamplerBinding("directLightSampler", directLightImageView, screenspaceSampler.get()));
        auto albedoImageView = passContext.GetImageView(this->viewportResources->albedo.imageViewProxy->Id());
        imageSamplerBindings.push_back(shaderDataSetInfo->MakeImageSamplerBinding("albedoSampler", albedoImageView, screenspaceSampler.get()));
        auto indirectLightView = passContext.GetImageView(this->viewportResources->indirectLight.imageViewProxy->Id());
        imageSamplerBindings.push_back(shaderDataSetInfo->MakeImageSamplerBinding("indirectLightSampler", indirectLightView, screenspaceSampler.get()));*/

        auto shaderDataSetBindings = legit::DescriptorSetBindings()
          .SetUniformBufferBindings(shaderData.uniformBufferBindings);
          //.SetImageSamplerBindings(imageSamplerBindings);
          //.SetStorageImageBindings(storageImageBindings);*/

        auto shaderDataSet = rendererData.core->GetDescriptorSetCache()->GetDescriptorSet(*shaderDataSetInfo, shaderDataSetBindings);

        passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, ShaderDataSetIndex, { shaderDataSet }, { shaderData.dynamicOffset });
        passContext.GetCommandBuffer().draw(4, 1, 0, 0);
      }
    }));
  }
}