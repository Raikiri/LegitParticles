#include "../Context/WindowData.h"

namespace almost
{

  WindowData InitWindow(legit::WindowFactory &windowFactory)
  {
    WindowData windowData;
    windowData.window = windowFactory.CreateWindow(1920, 1080, "Legit Vulkan renderer", nullptr, nullptr);

    return windowData;
  }
}