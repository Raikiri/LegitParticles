#pragma once
#include "../../Utils/GlmInclude.h"
#include <entity\registry.hpp>

namespace almost
{
  struct LinkComponent
  {
    entt::entity entities[2];
    float defLength;
  };
  struct LinkIndexComponent
  {
    size_t indices[2];
  };
}