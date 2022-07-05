#include "../Context/RendererData.h"
#include "../Context/ShaderTestbedRendererData.h"
#include "../Context/CameraData.h"
#include "lodepng.h"

#include <memory>
namespace almost
{

  void ReloadShader(almost::RendererData& rendererData, ShaderTestbedRendererData& shaderTestbedRendererData)
  {
    auto* core = rendererData.core.get();
    vk::PhysicalDevice physicalDevice = core->GetPhysicalDevice();
    vk::Device logicalDevice = core->GetLogicalDevice();

    shaderTestbedRendererData.shader.vertex.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/screenspaceQuad.vert.spv"));
    shaderTestbedRendererData.shader.fragment.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/animationTest.frag.spv"));
    shaderTestbedRendererData.shader.program.reset(new legit::ShaderProgram(shaderTestbedRendererData.shader.vertex.get(), shaderTestbedRendererData.shader.fragment.get()));
  }

  void ReloadImage(almost::RendererData& rendererData, ShaderTestbedRendererData& shaderTestbedRendererData)
  {
    auto* core = rendererData.core.get();
    shaderTestbedRendererData.sampler.reset(new legit::Sampler(core->GetLogicalDevice(), vk::SamplerAddressMode::eClampToEdge, vk::Filter::eNearest, vk::SamplerMipmapMode::eLinear));
    
    {
      std::vector<unsigned char> srcData;
      glm::uvec2 srcSize;
      auto res = lodepng::decode(srcData, srcSize.x, srcSize.y, "../data/Out/atlas.png");


      legit::ImageTexelData texelData = legit::CreateSimpleImageTexelData((glm::uint8*)srcData.data(), srcSize.x, srcSize.y, vk::Format::eR8G8B8A8Unorm);

      //auto texelData = legit::LoadKtxFromFile("../data/Textures/Cubemaps/5_16f.ktx");
      auto imageCreateDesc = legit::Image::CreateInfo2d(texelData.baseSize, uint32_t(texelData.mips.size()), uint32_t(texelData.layersCount), texelData.format, vk::ImageUsageFlagBits::eSampled | vk::ImageUsageFlagBits::eTransferDst);
      shaderTestbedRendererData.image = std::unique_ptr<legit::Image>(new legit::Image(core->GetPhysicalDevice(), core->GetLogicalDevice(), imageCreateDesc));
      legit::LoadTexelData(core, &texelData, shaderTestbedRendererData.image->GetImageData());
      shaderTestbedRendererData.imageView = std::unique_ptr<legit::ImageView>(
        new legit::ImageView(core->GetLogicalDevice(), shaderTestbedRendererData.image->GetImageData(), 0, 1, 0, 1));
    }
  }
  ShaderTestbedRendererData InitShaderTestbedRendererData(almost::RendererData& rendererData)
  {

    ShaderTestbedRendererData shaderTestbedRendererData;
    ReloadShader(rendererData, shaderTestbedRendererData);
    ReloadImage(rendererData, shaderTestbedRendererData);

    return shaderTestbedRendererData;
  }

  void RenderShaderTestbed(const almost::WindowData& windowData, almost::ShaderTestbedRendererData& testbedRendererData, almost::RendererData& rendererData, const almost::CameraData &cameraData)
  {
    static float val = 0.0f;
    ImGui::SliderFloat("Shader value", &val, 0.0f, 1.0f);
    #pragma pack(push, 1)
    struct GlobalDataBuffer
    {
      glm::vec2 camPos;
      glm::vec2 camFov;
      glm::vec2 viewportSize;
      float camAngle;
      float debugVal;
    };
    #pragma pack(pop)
    if (glfwGetKey(windowData.window, GLFW_KEY_V))
    {
      ReloadShader(rendererData, testbedRendererData);
    }
    rendererData.core->GetRenderGraph()->AddPass(legit::RenderGraph::RenderPassDesc()
      .SetColorAttachments({ { rendererData.frameInfo.swapchainImageViewProxyId, vk::AttachmentLoadOp::eClear, rendererData.clearColor}  })
      .SetProfilerInfo(legit::Colors::turqoise, "Grid")
      .SetRenderAreaExtent(rendererData.inFlightQueue->GetImageSize())
      .SetRecordFunc( [&testbedRendererData, &rendererData, &cameraData](legit::RenderGraph::RenderPassContext passContext)
    {
      const static uint32_t ShaderDataSetIndex = 0;
      const static uint32_t DrawCallDataSetIndex = 1;

      auto shader = testbedRendererData.shader.program.get();

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
          shaderDataBuffer->debugVal = val;
        }
        rendererData.frameInfo.memoryPool->EndSet();

        std::vector<legit::ImageSamplerBinding> imageSamplerBindings;
        imageSamplerBindings.push_back(shaderDataSetInfo->MakeImageSamplerBinding("animationTex", testbedRendererData.imageView.get(), testbedRendererData.sampler.get()));

        /*auto albedoImageView = passContext.GetImageView(this->viewportResources->albedo.imageViewProxy->Id());
        imageSamplerBindings.push_back(shaderDataSetInfo->MakeImageSamplerBinding("albedoSampler", albedoImageView, screenspaceSampler.get()));*/

        auto shaderDataSetBindings = legit::DescriptorSetBindings()
          .SetUniformBufferBindings(shaderData.uniformBufferBindings)
          .SetImageSamplerBindings(imageSamplerBindings);
          //.SetStorageImageBindings(storageImageBindings);*/

        auto shaderDataSet = rendererData.core->GetDescriptorSetCache()->GetDescriptorSet(*shaderDataSetInfo, shaderDataSetBindings);

        passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, ShaderDataSetIndex, { shaderDataSet }, { shaderData.dynamicOffset });
        passContext.GetCommandBuffer().draw(4, 1, 0, 0);
      }
    }));
  }
}