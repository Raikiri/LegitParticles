#pragma once
#include <memory>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#define NOMINMAX

#if defined(WIN32)
  #define GLFW_EXPOSE_NATIVE_WIN32
#else
  #define GLFW_EXPOSE_NATIVE_WAYLAND
#endif

#include <GLFW/glfw3native.h>
#include "../LegitVulkan/WindowDesc.h"
#include <iostream>
#include <string>

namespace legit
{
  struct WindowFactory
  {
    WindowFactory()
    {
      glfwSetErrorCallback([](int err_code, const char *desc)
      {
        std::cout << "GLFW error: " << desc << "\n";
      });
      #if !defined(WIN32)
        glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_WAYLAND);
      #endif
      if(!glfwInit())
      {
        std::cout << "GLFW failed to init!\n";
      }
      glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    }
    ~WindowFactory()
    {
      glfwTerminate();
    }
    
    struct Window
    {
    private:
      Window(size_t width, size_t height, std::string caption, GLFWmonitor *mon, GLFWwindow *share)
      {
        this->glfw_window = glfwCreateWindow(width, height, caption.c_str(), mon, nullptr);
      }
    public:
      ~Window()
      {
        glfwDestroyWindow(glfw_window);
      }

      legit::WindowDesc GetWindowDesc()
      {
        legit::WindowDesc windowDesc = {};
        #if defined(WIN32)
          windowDesc.hInstance = GetModuleHandle(NULL);
          windowDesc.hWnd = glfwGetWin32Window(window);
        #else
          windowDesc.display = glfwGetWaylandDisplay();
          windowDesc.surface = glfwGetWaylandWindow(glfw_window);
        #endif
        return windowDesc;
      }

      GLFWwindow* glfw_window;
      friend struct WindowFactory; 
    };
    
    std::unique_ptr<Window> CreateWindow(size_t width, size_t height, std::string caption, GLFWmonitor *mon = nullptr, GLFWwindow *share = nullptr)
    {
      return std::unique_ptr<Window>(new Window(width, height, caption, mon, share));
    }
  };
}