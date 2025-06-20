#pragma once

class DebugRenderer
{
public:
  DebugRenderer(legit::Core *_core)
  {
    this->core = _core;
    imageSpaceSampler.reset(new legit::Sampler(core->GetLogicalDevice(), vk::SamplerAddressMode::eClampToEdge, vk::Filter::eLinear, vk::SamplerMipmapMode::eNearest));
    nearestSampler.reset(new legit::Sampler(core->GetLogicalDevice(), vk::SamplerAddressMode::eClampToEdge, vk::Filter::eNearest, vk::SamplerMipmapMode::eNearest));
    ReloadShaders();
  }

  void RenderImageViews(legit::RenderGraph* renderGraph, legit::ShaderMemoryPool* memoryPool, legit::RenderGraph::ImageViewProxyId targetProxyId, std::vector<legit::RenderGraph::ImageViewProxyId> debugProxies, std::vector<legit::RenderGraph::ImageViewProxyId> debugUProxies = {})
  {
    glm::uvec2 viewportSize = renderGraph->GetMipSize(targetProxyId, 0);
    vk::Extent2D viewportExtent(viewportSize.x, viewportSize.y);

    glm::vec2 tileSize(0.3f, 0.3f);
    glm::vec2 tilePadding(0.02f, 0.02f);

    renderGraph->AddPass(legit::RenderGraph::RenderPassDesc()
      .SetColorAttachments({ {targetProxyId , vk::AttachmentLoadOp::eLoad} })
      .SetInputImages(std::vector<legit::RenderGraph::ImageViewProxyId>(debugProxies))
      .SetRenderAreaExtent(viewportExtent)
      .SetProfilerInfo(legit::Colors::clouds, "DebugInfoPass")
      .SetRecordFunc([this, tileSize, tilePadding, memoryPool, debugProxies, viewportSize](legit::RenderGraph::RenderPassContext passContext)
    {
      const auto& shader = debugRendererShader;
      auto pipeineInfo = this->core->GetPipelineCache()->BindGraphicsPipeline(passContext.GetCommandBuffer(), passContext.GetRenderPass()->GetHandle(), legit::DepthSettings::Disabled(), { legit::BlendSettings::Opaque() }, legit::VertexDeclaration(), vk::PrimitiveTopology::eTriangleFan, shader.program.get());
      {
        glm::vec2 currMin = tilePadding;

        for (auto debugProxyId : debugProxies)
        {
          const legit::DescriptorSetLayoutKey *vertexDataSetInfo = shader.program->GetSetInfo(VertexDataSetIndex);
          auto vertexData = memoryPool->BeginSet(vertexDataSetInfo);
          {
            auto quadDataBuffer = memoryPool->GetUniformBufferData<DebugRendererShader::QuadData>("QuadData");
            quadDataBuffer->minmax = glm::vec4(currMin.x, currMin.y, currMin.x + tileSize.x, currMin.y + tileSize.y);
          }
          memoryPool->EndSet();
          auto vertexDataSet = this->core->GetDescriptorSetCache()->GetDescriptorSet(*vertexDataSetInfo, vertexData.uniformBufferBindings, {}, {});

          const legit::DescriptorSetLayoutKey *fragmentDataSetInfo = shader.program->GetSetInfo(FragmentDataSetIndex);
          std::vector<legit::ImageSamplerBinding> imageSamplerBindings;
          imageSamplerBindings.push_back(fragmentDataSetInfo->MakeImageSamplerBinding("srcSampler", passContext.GetImageView(debugProxyId), imageSpaceSampler.get()));

          auto fragmentDataSet = this->core->GetDescriptorSetCache()->GetDescriptorSet(*fragmentDataSetInfo, {}, {}, { imageSamplerBindings });

          //passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, VertexDataSetIndex, { vertexDataSet, fragmentDataSet }, { vertexDataOffset });
          passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, VertexDataSetIndex, { vertexDataSet }, { vertexData.dynamicOffset });
          passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, FragmentDataSetIndex, { fragmentDataSet }, {});

          passContext.GetCommandBuffer().draw(4, 1, 0, 0);

          currMin.x += tileSize.x + tilePadding.x;
          if (currMin.x + tileSize.x > 1.0f)
          {
            currMin.x = tilePadding.x;
            currMin.y += tileSize.y + tilePadding.y;
          }
        }
      }
    }));

    size_t skipCount = debugProxies.size();
    renderGraph->AddPass(legit::RenderGraph::RenderPassDesc()
      .SetColorAttachments({ {targetProxyId , vk::AttachmentLoadOp::eLoad} })
      .SetInputImages(std::vector<legit::RenderGraph::ImageViewProxyId>(debugUProxies))
      .SetRenderAreaExtent(viewportExtent)
      .SetProfilerInfo(legit::Colors::clouds, "DebugInfoPass")
      .SetRecordFunc([this, skipCount, tileSize, tilePadding, memoryPool, debugUProxies, viewportSize](legit::RenderGraph::RenderPassContext passContext)
    {
      const auto& shader = debugURendererShader;
      auto pipeineInfo = this->core->GetPipelineCache()->BindGraphicsPipeline(passContext.GetCommandBuffer(), passContext.GetRenderPass()->GetHandle(), legit::DepthSettings::Disabled(), { legit::BlendSettings::Opaque() }, legit::VertexDeclaration(), vk::PrimitiveTopology::eTriangleFan, shader.program.get());
      {
        glm::vec2 currMin = tilePadding;
        for (size_t i = 0; i < skipCount; i++)
        {
          currMin.x += tileSize.x + tilePadding.x;
          if (currMin.x + tileSize.x > 1.0f)
          {
            currMin.x = tilePadding.x;
            currMin.y += tileSize.y + tilePadding.y;
          }
        }

        for (auto debugProxyId : debugUProxies)
        {
          const legit::DescriptorSetLayoutKey *vertexDataSetInfo = shader.program->GetSetInfo(VertexDataSetIndex);
          auto vertexData = memoryPool->BeginSet(vertexDataSetInfo);
          {
            auto quadDataBuffer = memoryPool->GetUniformBufferData<DebugRendererShader::QuadData>("QuadData");
            quadDataBuffer->minmax = glm::vec4(currMin.x, currMin.y, currMin.x + tileSize.x, currMin.y + tileSize.y);
          }
          memoryPool->EndSet();
          auto vertexDataSet = this->core->GetDescriptorSetCache()->GetDescriptorSet(*vertexDataSetInfo, vertexData.uniformBufferBindings, {}, {});

          const legit::DescriptorSetLayoutKey *fragmentDataSetInfo = shader.program->GetSetInfo(FragmentDataSetIndex);
          std::vector<legit::ImageSamplerBinding> imageSamplerBindings;
          imageSamplerBindings.push_back(fragmentDataSetInfo->MakeImageSamplerBinding("srcUSampler", passContext.GetImageView(debugProxyId), nearestSampler.get()));

          auto fragmentDataSet = this->core->GetDescriptorSetCache()->GetDescriptorSet(*fragmentDataSetInfo, {}, {}, { imageSamplerBindings });

          //passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, VertexDataSetIndex, { vertexDataSet, fragmentDataSet }, { vertexDataOffset });
          passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, VertexDataSetIndex, { vertexDataSet }, { vertexData.dynamicOffset });
          passContext.GetCommandBuffer().bindDescriptorSets(vk::PipelineBindPoint::eGraphics, pipeineInfo.pipelineLayout, FragmentDataSetIndex, { fragmentDataSet }, {});

          passContext.GetCommandBuffer().draw(4, 1, 0, 0);

          currMin.x += tileSize.x + tilePadding.x;
          if (currMin.x + tileSize.x > 1.0f)
          {
            currMin.x = tilePadding.x;
            currMin.y += tileSize.y + tilePadding.y;
          }
        }
      }
    }));
  }

  void ReloadShaders()
  {
    debugRendererShader.vertex.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/Common/debugRenderer.vert.spv"));
    debugRendererShader.fragment.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/Common/debugRenderer.frag.spv"));
    debugRendererShader.program.reset(new legit::ShaderProgram(debugRendererShader.vertex.get(), debugRendererShader.fragment.get()));

    debugURendererShader.vertex.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/Common/debugRenderer.vert.spv"));
    debugURendererShader.fragment.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/Common/debugURenderer.frag.spv"));
    debugURendererShader.program.reset(new legit::ShaderProgram(debugURendererShader.vertex.get(), debugURendererShader.fragment.get()));
  }
private:


  const static uint32_t VertexDataSetIndex = 0;
  const static uint32_t FragmentDataSetIndex = 1;

  struct DebugRendererShader
  {
#pragma pack(push, 1)
    struct QuadData
    {
      glm::vec4 minmax;
    };
#pragma pack(pop)

    std::unique_ptr<legit::Shader> vertex;
    std::unique_ptr<legit::Shader> fragment;
    std::unique_ptr<legit::ShaderProgram> program;
  } debugRendererShader;

  struct DebugURendererShader
  {
#pragma pack(push, 1)
    struct QuadData
    {
      glm::vec4 minmax;
    };
#pragma pack(pop)

    std::unique_ptr<legit::Shader> vertex;
    std::unique_ptr<legit::Shader> fragment;
    std::unique_ptr<legit::ShaderProgram> program;
  } debugURendererShader;

  vk::Extent2D viewportSize;

  std::unique_ptr<legit::Sampler> imageSpaceSampler;
  std::unique_ptr<legit::Sampler> nearestSampler;

  legit::Core *core;
};