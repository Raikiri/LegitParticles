#pragma once
#include <entity/registry.hpp>
#include "../../Utils/GlmInclude.h"


#define POSITION_BASED
namespace almost
{
  struct ParticleComponent
  {
    glm::vec2 pos;
    #if defined(POSITION_BASED)
      glm::vec2 prevPos;
    #else
      glm::vec2 velocity;
    #endif
    glm::vec2 acceleration;
  };

  struct Influence
  {
    entt::entity particleEntity;
    float weight;
  };
  struct CoarseMultigridComponent
  {
    std::vector<Influence> influences;
  };

  struct FineMultigridComponent
  {
    std::vector<Influence> influences;
  };
  struct ParticleIndexComponent
  {
    size_t index;
  };
  struct MassComponent
  {
    float invMass;
    bool usesGravity;
  };
  struct DefPosComponent
  {
    glm::vec2 defPos;
    bool isDraggable;
  };
}