#include "../Context/RendererData.h"
#include "../Context/MeshRendererData.h"
#include "../Context/CameraData.h"

#include <memory>
namespace almost
{

  MeshRendererData InitMeshRendererData(almost::RendererData& rendererData)
  {
    auto *core = rendererData.core.get();
    vk::PhysicalDevice physicalDevice = core->GetPhysicalDevice();
    vk::Device logicalDevice = core->GetLogicalDevice();

    MeshRendererData meshRendererData;
    meshRendererData.maxVerticesCount = 1000000;
    meshRendererData.maxIndicesCount = 1000000;
    meshRendererData.vertexData = nullptr;
    meshRendererData.indexData = nullptr;

    meshRendererData.vertexBuffer = std::unique_ptr<legit::Buffer>(new legit::Buffer(core->GetPhysicalDevice(), core->GetLogicalDevice(), sizeof(MeshRendererData::Vertex) * meshRendererData.maxVerticesCount, vk::BufferUsageFlagBits::eVertexBuffer, vk::MemoryPropertyFlagBits::eHostCoherent));
    meshRendererData.indexBuffer = std::unique_ptr<legit::Buffer>(new legit::Buffer(core->GetPhysicalDevice(), core->GetLogicalDevice(), sizeof(MeshRendererData::Index) * meshRendererData.maxIndicesCount, vk::BufferUsageFlagBits::eIndexBuffer, vk::MemoryPropertyFlagBits::eHostCoherent));

    meshRendererData.shader.vertex.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/meshRenderer.vert.spv"));
    meshRendererData.shader.fragment.reset(new legit::Shader(core->GetLogicalDevice(), "../data/Shaders/spirv/meshRenderer.frag.spv"));
    meshRendererData.shader.program.reset(new legit::ShaderProgram(meshRendererData.shader.vertex.get(), meshRendererData.shader.fragment.get()));


    legit::VertexDeclaration vertexDecl;
    vertexDecl.AddVertexInputBinding(0, sizeof(MeshRendererData::Vertex));
    vertexDecl.AddVertexAttribute(0, offsetof(MeshRendererData::Vertex, pos), legit::VertexDeclaration::AttribTypes::vec3, 0);
    vertexDecl.AddVertexAttribute(0, offsetof(MeshRendererData::Vertex, uv), legit::VertexDeclaration::AttribTypes::vec2, 1);
    vertexDecl.AddVertexAttribute(0, offsetof(MeshRendererData::Vertex, color), legit::VertexDeclaration::AttribTypes::vec4, 2);
    meshRendererData.vertexDecl = vertexDecl;

    return meshRendererData;
  }

  void MapMeshRenderer(almost::MeshRendererData& meshRendererData)
  {
    meshRendererData.vertexData = (MeshRendererData::Vertex*)meshRendererData.vertexBuffer->Map();
    meshRendererData.indexData = (MeshRendererData::Index*)meshRendererData.indexBuffer->Map();
  }
  void UnmapMeshRenderer(almost::MeshRendererData& meshRendererData)
  {
    meshRendererData.vertexBuffer->Unmap();
    meshRendererData.indexBuffer->Unmap();
  }


  void RenderMeshes(const almost::WindowData& windowData, const almost::MeshRendererData& meshRendererData, almost::RendererData& rendererData, const almost::CameraData &cameraData)
  {
    #pragma pack(push, 1)
    struct GlobalDataBuffer
    {
      glm::mat4 viewProjMatrix;
    };
    #pragma pack(pop)
    glm::mat4 viewProjMatrix = almost::GetViewProjMatrix(cameraData.pos, cameraData.ang, cameraData.fov);

    rendererData.core->GetRenderGraph()->AddPass(legit::RenderGraph::RenderPassDesc()
      .SetColorAttachments({ { rendererData.frameInfo.swapchainImageViewProxyId, vk::AttachmentLoadOp::eLoad}  })
      .SetProfilerInfo(legit::Colors::pomegranate, "Meshes")
      .SetRenderAreaExtent(rendererData.inFlightQueue->GetImageSize())
      .SetRecordFunc( [viewProjMatrix, &meshRendererData, &rendererData](legit::RenderGraph::RenderPassContext passContext)
    {
      const static uint32_t ShaderDataSetIndex = 0;
      const static uint32_t DrawCallDataSetIndex = 1;

      auto shader = meshRendererData.shader.program.get();

      auto pipeineInfo = rendererData.core->GetPipelineCache()->BindGraphicsPipeline(passContext.GetCommandBuffer(), passContext.GetRenderPass()->GetHandle(), legit::DepthSettings::Disabled(), { legit::BlendSettings::AlphaBlend() }, meshRendererData.vertexDecl, vk::PrimitiveTopology::eTriangleList, shader);
      {
        const legit::DescriptorSetLayoutKey *shaderDataSetInfo = shader->GetSetInfo(ShaderDataSetIndex);
        auto shaderData = rendererData.frameInfo.memoryPool->BeginSet(shaderDataSetInfo);
        {
          auto shaderDataBuffer = rendererData.frameInfo.memoryPool->GetUniformBufferData<GlobalDataBuffer>("GlobalDataBuffer");

          shaderDataBuffer->viewProjMatrix = viewProjMatrix;
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
        passContext.GetCommandBuffer().bindVertexBuffers(0, { meshRendererData.vertexBuffer->GetHandle() }, { 0 });
        passContext.GetCommandBuffer().bindIndexBuffer(meshRendererData.indexBuffer->GetHandle(), 0, vk::IndexType::eUint32);
        passContext.GetCommandBuffer().drawIndexed(glm::uint32_t(meshRendererData.indicesCount), 1, 0, 0, 0);
      }
    }));
  }
}