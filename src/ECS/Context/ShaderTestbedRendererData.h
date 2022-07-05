#pragma once
#include <vector>
#include "../../Utils/GlmInclude.h"
#include "../../LegitVulkan/LegitVulkan.h"

namespace almost
{
  struct ShaderTestbedRendererData
  {
    struct Shader
    {
      std::unique_ptr<legit::Shader> vertex;
      std::unique_ptr<legit::Shader> fragment;
      std::unique_ptr<legit::ShaderProgram> program;

    } shader;

    std::unique_ptr<legit::Sampler> sampler;

    std::unique_ptr<legit::Image> image;
    std::unique_ptr<legit::ImageView> imageView;
  };
}