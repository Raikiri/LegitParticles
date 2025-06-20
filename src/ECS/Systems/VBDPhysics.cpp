#include <string>
#include "PhysicsCommon.h"
#include "VBDPhysics.h"
#include "../Context/MeshRendererData.h"
#include "../Context/InputData.h"
#include "../Context/CameraData.h"
#include "../Context/WindowData.h"
#include "../../Maths/StrainConstraints/StrainConstraints.h"
#include "imgui.h"
//#include "../../Maths/Tensor2.h"

namespace almost
{
  void ProcessPhysicsVBD(
    std::vector<entt::registry>& regLayers,
    legit::CpuProfiler& profiler)
  {
    auto &reg = regLayers[0];
    float dt = 1e-2f;
    {
      auto physicsTask = profiler.StartScopedTask("[Physics] PreStep", legit::Colors::sunFlower);
      PreStep(
        almost::ParticleGroup::Get(reg),
        almost::LinkGroup::Get(reg),
        almost::TriangleGroup::Get(reg));
    }
    static bool use_simple_solver = true;
    ImGui::Checkbox("Use simple solver", &use_simple_solver);
    
    if(use_simple_solver)
    {
      {
        auto physicsTask = profiler.StartScopedTask("[Physics] Links", legit::Colors::orange);
        auto particleGroup = almost::ParticleGroup::Get(reg);
        auto linkGroup = almost::LinkGroup::Get(reg);
        for(int i = 0; i < 10; i++)
        {
          ProjectLinks(
            particleGroup.raw<ParticleComponent>(), particleGroup.raw<MassComponent>(), particleGroup.size(),
            linkGroup.raw<LinkComponent>(), linkGroup.raw<LinkIndexComponent>(), linkGroup.size());
        }
      }
      {
        auto physicsTask = profiler.StartScopedTask("[Physics] Integration", legit::Colors::orange);

        auto particleGroup = almost::ParticleGroup::Get(reg);
        IntegrateParticles(particleGroup.raw<ParticleComponent>(), particleGroup.size(), dt);
      }
    }
  }
}