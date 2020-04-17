#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../LegitVulkan/LegitVulkan.h"

namespace almost
{
  struct FullscreenRendererData
  {
    struct Shader
    {
      std::unique_ptr<legit::Shader> vertex;
      std::unique_ptr<legit::Shader> fragment;
      std::unique_ptr<legit::ShaderProgram> program;
    } shader;
  };
}