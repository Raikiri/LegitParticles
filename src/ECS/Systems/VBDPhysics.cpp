#include <string>
#include "VBDPhysics.h"
#include "../Context/MeshRendererData.h"
#include "../Context/InputData.h"
#include "../Context/CameraData.h"
#include "../Context/WindowData.h"
#include "../../Maths/StrainConstraints/StrainConstraints.h"
//#include "../../Maths/Tensor2.h"

namespace almost
{
  void ProcessPhysicsVBD(
    std::vector<entt::registry>& regLayers,
    legit::CpuProfiler& profiler)
  {
    auto &reg = regLayers[0];
    float dt = 1e-2f;
    PreStep(
      almost::ParticleGroup::Get(reg),
      almost::LinkGroup::Get(reg),
      almost::TriangleGroup::Get(reg),
      profiler);

    auto particleGroup0 = almost::ParticleGroup::Get(reg);
    IntegrateParticles(particleGroup0.raw<ParticleComponent>(), particleGroup0.size(), dt);
  }
}