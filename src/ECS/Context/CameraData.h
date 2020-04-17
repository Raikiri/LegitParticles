#pragma once
#include "../../Utils/GlmInclude.h"

namespace almost
{
  struct CameraData
  {
    glm::vec2 pos;
    glm::vec2 fov;
    float ang;
  };

  glm::vec2 GetFov2(float horFov, glm::vec2 renderArea)
  {
    return glm::vec2(horFov, horFov / renderArea.x * renderArea.y);
  }
  glm::mat4 GetViewProjMatrix(glm::vec2 pos, float ang, glm::vec2 fov)
  {
    return glm::scale(glm::vec3(1.0f, -1.0f, 1.0f)) * glm::scale(glm::vec3(1.0 / fov.x, 1.0f / fov.y, 1.0f)) * glm::rotate(-ang, glm::vec3(0.0f, 0.0f, 1.0f)) * glm::translate(glm::vec3(-pos, 0.0f));
  }

  glm::vec2 Project(glm::vec2 pos, float ang, glm::vec2 fov, glm::vec2 worldPos)
  {
    glm::vec2 rightVec = glm::vec2(cos(ang), sin(ang));
    glm::vec2 upVec = glm::vec2(-rightVec.y, rightVec.x);
    return glm::vec2(glm::dot(worldPos - pos, rightVec), glm::dot(worldPos - pos, upVec)) / glm::vec2(fov * 2.0f) + glm::vec2(0.5f);
  }

  glm::vec2 Unproject(glm::vec2 pos, float ang, glm::vec2 fov, glm::vec2 screenPos)
  {
    glm::vec2 rightVec = glm::vec2(cos(ang), sin(ang));
    glm::vec2 upVec = glm::vec2(-rightVec.y, rightVec.x);
    return pos + rightVec * fov.x * (screenPos.x - 0.5f) * 2.0f + upVec * fov.y * (screenPos.y - 0.5f) * 2.0f;
  }

  struct AABB2
  {
    glm::vec2 min;
    glm::vec2 max;
  };
  AABB2 GetCameraAABB(glm::vec2 pos, float ang, glm::vec2 fov)
  {
    AABB2 aabb;
    glm::vec2 rightVec = glm::vec2(cos(ang), sin(ang));
    glm::vec2 upVec = glm::vec2(-rightVec.y, rightVec.x);
    aabb.min = pos - rightVec * fov.x - upVec * fov.y;
    aabb.max = pos + rightVec * fov.x + upVec * fov.y;
    return aabb;
  }

}