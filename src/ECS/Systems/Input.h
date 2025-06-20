#include "../Context/WindowData.h"
#include "../Context/InputData.h"
namespace almost
{

  InputData InitInputData(const almost::WindowData& windowData)
  {
    InputData inputData;
    inputData.prevFrameTime = std::chrono::system_clock::now();
    inputData.deltaTime = 0.0f;
    glfwGetCursorPos(windowData.window->glfw_window, &inputData.mousePos.x, &inputData.mousePos.y);
    inputData.prevMousePos = inputData.mousePos;
    inputData.worldMousePos = glm::vec2(0.0f);
    inputData.prevWorldMousePos = glm::vec2(0.0f);
    inputData.isWindowClosed = false;
    return inputData;
  }
  void ProcessInput(const almost::WindowData& windowData, almost::InputData& inputData)
  {
    auto currFrameTime = std::chrono::system_clock::now();
    inputData.deltaTime = std::chrono::duration<float>(currFrameTime - inputData.prevFrameTime).count();
    inputData.prevFrameTime = currFrameTime;

    glfwPollEvents();
    inputData.prevMousePos = inputData.mousePos;
    glfwGetCursorPos(windowData.window->glfw_window, &inputData.mousePos.x, &inputData.mousePos.y);
    inputData.isWindowClosed = glfwWindowShouldClose(windowData.window->glfw_window);
  }

}