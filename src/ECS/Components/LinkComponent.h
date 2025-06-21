#pragma once
#include "../../Utils/GlmInclude.h"
#include <entity/registry.hpp>

namespace almost
{
  struct LinkComponent
  {
    entt::entity entities[2];
    float defLength = 0.0f;
    float stiffness = 0.0f;
    float lambda = 0.0f;
    float deltaC = 0.0f;
  };
  struct LinkIndexComponent
  {
    size_t GetOtherParticleIdx(size_t particle_idx) const
    {
      if(indices[0] == particle_idx)
      {
        return indices[1];
      }else
      {
        assert(indices[1] == particle_idx);
        return indices[0];
      }
    }

    size_t indices[2];
  };
}