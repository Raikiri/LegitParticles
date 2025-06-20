#pragma once

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#if defined(WIN32)
  #define GLFW_EXPOSE_NATIVE_WIN32
#else
  #define GLFW_EXPOSE_NATIVE_WAYLAND
#endif
#include <GLFW/glfw3native.h>

namespace almost
{
  struct WindowData
  {
    GLFWwindow* window;
  };
}