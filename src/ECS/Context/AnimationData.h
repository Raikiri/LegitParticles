#pragma once
#include <vector>
#include <string>
#include "../../Utils/GlmInclude.h"
#include "../../Maths/Tensor2.h"
#include "../../Maths/SpriterAnimations/SpriterImporter.h"
namespace almost
{
  struct AnimationData
  {
    spriter::Object spriterObject;
    size_t currKeyIndex = 0;

    FlattenedAnimation flattenedAnimation;
  };
}