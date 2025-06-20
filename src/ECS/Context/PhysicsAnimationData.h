#pragma once
#include <vector>
#include <string>
#include "../../Utils/GlmInclude.h"

namespace almost
{
  struct SmoothState
  {
    glm::vec2 position;
    glm::vec2 velocity;
  };

  struct PhysicsAnimationData
  {
    glm::vec2 objPos0;
    glm::vec2 prevObjPos0;
    glm::vec2 prevMouseWorldPos;
    bool isReset = false;

    std::chrono::system_clock::time_point lastUpdateTime;

    SmoothState objState1;
  };
}