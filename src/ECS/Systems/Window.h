#include "../Context/WindowData.h"

namespace almost
{

  WindowData InitWindow()
  {
    WindowData windowData;
    glfwInit();
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    windowData.window = glfwCreateWindow(1600, 1024, "Legit Vulkan renderer", nullptr, nullptr);

    return windowData;
  }
}