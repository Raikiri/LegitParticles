#pragma once
#include "../../LegitProfiler/ProfilerTask.h"
#include "../../LegitVulkan/CpuProfiler.h"

#include "../Context/PhysicsData.h"
//#include "../Context/MeshRendererData.h"
#include <entity/registry.hpp>
#include "../../Utils/GroupArg.h"
#include "../../Maths/StackStorage.h"
#include "PhysicsCommon.h"

namespace almost
{
  void ProcessPhysicsVBD(
    std::vector<entt::registry> &regLayers,
    legit::CpuProfiler& profiler);
}