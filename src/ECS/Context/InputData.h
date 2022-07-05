#pragma once
#include "../../Utils/GlmInclude.h"
#include <chrono>
#include <ctime>

namespace almost
{
  struct InputData
  {
    glm::f64vec2 mousePos;
    glm::f64vec2 prevMousePos;
    glm::vec2 worldMousePos;
    glm::vec2 prevWorldMousePos;
    float deltaTime;
    std::chrono::system_clock::time_point prevFrameTime;
    bool isWindowClosed;
  };
}