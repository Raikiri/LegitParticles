#pragma once
#include "../../LegitVulkan/LegitVulkan.h"
#include "../../Utils/GlmInclude.h"
#include "imgui.h"
#include "../../LegitProfiler/ImGuiProfilerRenderer.h"
#include "../../LegitImGui/ImGuiRenderer.h"

namespace almost
{
  struct RendererData
  {
    std::unique_ptr<legit::Core> core;
    std::unique_ptr<ImGuiRenderer> imguiRenderer;
    std::unique_ptr<legit::InFlightQueue> inFlightQueue;
    legit::WindowDesc windowDesc;
    vk::ClearColorValue clearColor;
    legit::InFlightQueue::FrameInfo frameInfo;
    ImGuiUtils::ProfilersWindow profilersWindow;
  };
}