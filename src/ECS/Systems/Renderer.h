#include "../Context/RendererData.h"
#include "../Context/WindowData.h"
#include "../Context/MeshRendererData.h"
namespace almost
{
  glm::uvec2 GetGlfwWindowClientSize(GLFWwindow *window)
  {
    int width = 0, height = 0;
    glfwGetWindowSize(window, &width, &height);
    return {width, height};
  }

  RendererData InitRenderer(almost::WindowData& windowData)
  {
    RendererData rendererData;
    //const char* glfwExtensions[] = { "VK_KHR_surface", "VK_KHR_win32_surface" };
    //uint32_t glfwExtensionCount = 2;

    rendererData.windowDesc = windowData.window->GetWindowDesc();

    uint32_t glfwExtensionCount = 0;
    const char** glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);
    rendererData.core = std::make_unique<legit::Core>(glfwExtensions, glfwExtensionCount, &rendererData.windowDesc, true);

    rendererData.imguiRenderer = std::make_unique<ImGuiRenderer>(rendererData.core.get(), windowData.window->glfw_window);
    float gray = 0.02f;
    rendererData.clearColor = vk::ClearColorValue(std::array<float, 4>{gray, gray, gray, 1.0f});

    return rendererData;
  }

  void StartFrame(const almost::WindowData& windowData, almost::RendererData& rendererData)
  {
    rendererData.imguiRenderer->ProcessInput(windowData.window->glfw_window);
    if (!rendererData.inFlightQueue)
    {
      std::cout << "recreated\n";
      rendererData.core->ClearCaches();
      //core->GetRenderGraph()->Clear();
      rendererData.inFlightQueue = std::unique_ptr<legit::InFlightQueue>(new legit::InFlightQueue(rendererData.core.get(), rendererData.windowDesc, GetGlfwWindowClientSize(windowData.window->glfw_window), 2, vk::PresentModeKHR::eMailbox));
      rendererData.imguiRenderer->RecreateSwapchainResources(rendererData.inFlightQueue->GetImageSize(), rendererData.inFlightQueue->GetInFlightFramesCount());

      auto& imguiIO = ImGui::GetIO();
      imguiIO.DeltaTime = 1.0f / 60.0f;              // set the time elapsed since the previous frame (in seconds)
      imguiIO.DisplaySize.x = float(rendererData.inFlightQueue->GetImageSize().width);             // set the current display width
      imguiIO.DisplaySize.y = float(rendererData.inFlightQueue->GetImageSize().height);             // set the current display height here
    }

    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_V))
    {
      //renderer.ReloadShaders();
    }
    rendererData.frameInfo = rendererData.inFlightQueue->BeginFrame();
    ImGui::NewFrame();
    auto& gpuProfilerData = rendererData.inFlightQueue->GetLastFrameGpuProfilerData();
    auto& cpuProfilerData = rendererData.inFlightQueue->GetLastFrameCpuProfilerData();

    //renderFunc(frameInfo);
    if (!rendererData.profilersWindow.stopProfiling)
    {
      auto profilersTask = rendererData.inFlightQueue->GetCpuProfiler().StartScopedTask("Prf processing", legit::Colors::sunFlower);

      rendererData.profilersWindow.gpuGraph.LoadFrameData(gpuProfilerData.data(), gpuProfilerData.size());
      rendererData.profilersWindow.cpuGraph.LoadFrameData(cpuProfilerData.data(), cpuProfilerData.size());
    }

    {
      auto profilersTask = rendererData.inFlightQueue->GetCpuProfiler().StartScopedTask("Prf rendering", legit::Colors::belizeHole);
      rendererData.profilersWindow.Render();
    }

    //ImGui::ShowStyleEditor();
    //ImGui::ShowDemoWindow();
  }
  void FinishFrame(const almost::WindowData& windowData, almost::RendererData& rendererData)
  {
    ImGui::EndFrame();

    ImGui::Render();
    rendererData.imguiRenderer->RenderFrame(rendererData.frameInfo, windowData.window->glfw_window, ImGui::GetDrawData());

    try
    {
      rendererData.inFlightQueue->EndFrame();
    }
    catch (vk::OutOfDateKHRError err)
    {
      rendererData.core->WaitIdle();
      rendererData.inFlightQueue.reset();
    }
    
    int width = 0, height = 0;
    glfwGetWindowSize(windowData.window->glfw_window, &width, &height);
    if(rendererData.inFlightQueue && rendererData.inFlightQueue->GetImageSize().width != width || rendererData.inFlightQueue->GetImageSize().height != height)
    {
      rendererData.core->WaitIdle();
      rendererData.inFlightQueue.reset();
    }

  }
}