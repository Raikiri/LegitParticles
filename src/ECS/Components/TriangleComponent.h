#pragma once
#include "../../Utils/GlmInclude.h"
#include <entity\registry.hpp>
#include "../../Maths/StrainConstraints/StrainConstraints.h"

namespace almost
{
  struct TriangleComponent
  {
    entt::entity entities[3];
    almost::StaticTensor<float, StrainDynamics::UvDim, StrainDynamics::RefDim> uvFromRef;
  };
  struct TriangleIndexComponent
  {
    size_t indices[3];
  };
}