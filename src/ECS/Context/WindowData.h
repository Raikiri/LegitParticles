#pragma once
#include "LegitGLFW/LegitGLFW.h"

namespace almost
{
  struct WindowData
  {
    std::unique_ptr<legit::WindowFactory::Window> window;
  };
}