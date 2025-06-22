#pragma once
#include <vector>

#include "../../Utils/GlmInclude.h"
#include "LegitVulkan/LegitVulkan.h"

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

  static void SubmitLine(almost::MeshRendererData& meshRendererData, glm::vec2 point0, glm::vec2 point1, float width, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    glm::vec2 dir = glm::normalize(point1 - point0);
    glm::vec2 tangent = glm::vec2(-dir.y, dir.x);

    meshRendererData.vertexData[vertexOffset + 0] = { glm::vec3(point0 + tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 1] = { glm::vec3(point0 - tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 2] = { glm::vec3(point1 - tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 3] = { glm::vec3(point1 + tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.indexData[indexOffset + 0] = glm::uint32_t(vertexOffset + 0);
    meshRendererData.indexData[indexOffset + 1] = glm::uint32_t(vertexOffset + 1);
    meshRendererData.indexData[indexOffset + 2] = glm::uint32_t(vertexOffset + 2);

    meshRendererData.indexData[indexOffset + 3] = glm::uint32_t(vertexOffset + 0);
    meshRendererData.indexData[indexOffset + 4] = glm::uint32_t(vertexOffset + 2);
    meshRendererData.indexData[indexOffset + 5] = glm::uint32_t(vertexOffset + 3);
    meshRendererData.verticesCount += 4;
    meshRendererData.indicesCount += 6;
  }

  static void SubmitCircle(almost::MeshRendererData& meshRendererData, glm::vec2 center, float radius, size_t sectorsCount, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    for (size_t sectorIndex = 0; sectorIndex < sectorsCount; sectorIndex++)
    {
      float ang = (float(sectorIndex) / float(sectorsCount)) * 2.0f * glm::pi<float>();
      meshRendererData.vertexData[vertexOffset + sectorIndex] = { glm::vec3(center + glm::vec2(cos(ang), sin(ang)) * radius, 0.0f), glm::vec2(0.0f, 0.0f), color };

      if (sectorIndex < sectorsCount - 2)
        meshRendererData.indexData[indexOffset + sectorIndex * 3 + 0] = glm::uint32_t(vertexOffset);
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 1] = glm::uint32_t(vertexOffset + sectorIndex + 1);
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 2] = glm::uint32_t(vertexOffset + sectorIndex + 2);
    }
    meshRendererData.verticesCount += sectorsCount;
    meshRendererData.indicesCount += (sectorsCount - 2) * 3;

  }

  static void SubmitTriangle(almost::MeshRendererData& meshRendererData, std::array<glm::vec2, 3> points, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    for (size_t vertexIndex = 0; vertexIndex < 3; vertexIndex++)
    {
      meshRendererData.vertexData[vertexOffset + vertexIndex] = { glm::vec3(points[vertexIndex], 0.0f), glm::vec2(0.0f, 0.0f), color };

      meshRendererData.indexData[indexOffset + vertexIndex] = glm::uint32_t(vertexOffset + vertexIndex);
    }
    meshRendererData.verticesCount += 3;
    meshRendererData.indicesCount += 3;
  }
}