#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../LegitVulkan/LegitVulkan.h"

namespace almost
{
  struct MeshRendererData
  {
    #pragma pack(push, 1)
    struct Vertex
    {
      glm::vec3 pos;
      glm::vec2 uv;
      glm::vec4 color;
    };
    #pragma pack(pop)
    using Index = glm::uint32_t;

    legit::VertexDeclaration vertexDecl;
    std::unique_ptr<legit::Buffer> vertexBuffer;
    Vertex* vertexData;
    std::unique_ptr<legit::Buffer> indexBuffer;
    Index* indexData;
    size_t indicesCount;
    size_t verticesCount;
    size_t maxVerticesCount;
    size_t maxIndicesCount;

    struct Shader
    {
      std::unique_ptr<legit::Shader> vertex;
      std::unique_ptr<legit::Shader> fragment;
      std::unique_ptr<legit::ShaderProgram> program;
    } shader;
  };

}