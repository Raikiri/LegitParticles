#include <string>
#include "PhysicsCommon.h"
#include "VBDPhysics.h"
#include "../Context/MeshRendererData.h"
#include "../Context/InputData.h"
#include "../Context/CameraData.h"
#include "../Context/WindowData.h"
#include "../../Maths/StrainConstraints/StrainConstraints.h"
#include "imgui.h"
#include "PhysicsGraphs.h"
//#include "../../Maths/Tensor2.h"

namespace almost
{
  void GetParticleInertialPositions(const ParticleComponent* particle_components, glm::vec2* inertial_positions, size_t particles_count, float dt)
  {
    for (size_t particle_idx = 0; particle_idx < particles_count; particle_idx++)
    {
      const ParticleComponent& particle_component = particle_components[particle_idx];
      #if defined(POSITION_BASED)
        inertial_positions[particle_idx] = particle_component.pos + (particle_component.pos - particle_component.prevPos) + particle_component.acceleration * dt * dt;
      #else
        inertial_positions[particle_idx] = particle_component.pos + particle_component.velocity * dt + particle_component.acceleration * dt * dt;
      #endif
    }
  }
  
  void AssignInertialPositions(ParticleComponent* particle_components, const glm::vec2* inertial_positions, size_t particles_count)
  {
    for (size_t particle_idx = 0; particle_idx < particles_count; particle_idx++)
    {
      ParticleComponent& particle_component = particle_components[particle_idx];
      particle_component.prevPos = particle_component.pos;
      particle_component.pos = inertial_positions[particle_idx];
    }
  }
  
  float GetLinkEnergy(
    glm::vec2 pos,
    size_t particle_idx,
    const ParticleComponent* particle_components,
    size_t particles_count,
    LinkIndexComponent link_particles,
    float def_len,
    float stiffness)
  {
    ConstraintGraph::ParticleIdx curr_particle_idx = particle_idx;
    ConstraintGraph::ParticleIdx other_particle_idx = ConstraintGraph::ParticleIdx(-1);
    if(link_particles.indices[0] == particle_idx)
    {
      other_particle_idx = link_particles.indices[1];
    }else
    {
      assert(link_particles.indices[1] == particle_idx);
      other_particle_idx = link_particles.indices[0];
    }
    float curr_len = glm::length(particle_components[other_particle_idx].pos - pos);
    float len_delta = curr_len - def_len;
    return len_delta * len_delta * stiffness * 0.5f;
  }
  
  float GetParticleEnergy(
    glm::vec2 pos,
    size_t particle_idx,
    float particle_mass,
    const ParticleComponent* particle_components,
    const glm::vec2 inertial_position,
    size_t particles_count,
    LinkIndexComponent *link_indices,
    LinkComponent *links,
    size_t links_count,
    const ConstraintGraph &constraint_graph,
    float dt)
  {
    float delta_len = glm::length(pos - inertial_position);
    float energy = delta_len * delta_len * 0.5f * particle_mass / (dt * dt);
    for(size_t link_num = 0; link_num < constraint_graph.particle_infos[particle_idx].links_count; link_num++)
    {
      ConstraintGraph::LinkIdx link_idx = constraint_graph.particle_links[constraint_graph.particle_infos[particle_idx].links_start + link_num];
      energy += GetLinkEnergy(
        pos,
        particle_idx,
        particle_components,
        particles_count,
        link_indices[link_idx],
        links[link_idx].defLength,
        1000.0f);
    }
    return energy;
  }
  void ProjectVBD(
    ParticleComponent* particle_components,
    const MassComponent* mass_components,
    const glm::vec2* inertial_positions,
    size_t particles_count,
    LinkIndexComponent *link_indices,
    LinkComponent *links,
    size_t links_count,
    const ConstraintGraph &constraint_graph,
    float dt)
  {
    for(size_t particle_idx = 0; particle_idx < particles_count; particle_idx++)
    {
      auto &particle = particle_components[particle_idx];
      auto inv_mass = mass_components[particle_idx].invMass;
      auto inertial_pos = inertial_positions[particle_idx];

      glm::vec2 best_pos = inertial_pos;
      float best_energy = 1e9f;
      glm::vec2 ratio;
      float step = 1e-2f;
      for(ratio.y = 0.0f; ratio.y < 1.0f; ratio.y += step)
      {
        for(ratio.x = 0.0f; ratio.x < 1.0f; ratio.x += step)
        {
          glm::vec2 test_pos = particle.pos + (ratio - glm::vec2(0.5f)) * 20.0f;
          float test_energy = GetParticleEnergy(
            test_pos,
            particle_idx,
            1.0f / std::max(1e-7f, inv_mass),
            particle_components,
            inertial_pos,
            particles_count,
            link_indices,
            links,
            links_count,
            constraint_graph,
            dt);
          if(test_energy < best_energy)
          {
            best_energy = test_energy;
            best_pos = test_pos;
          }
        }
      }
      particle.pos = best_pos;
    }
  }

  void ProcessPhysicsVBD(
    std::vector<entt::registry>& regLayers,
    legit::CpuProfiler& profiler)
  {
    auto &reg = regLayers[0];
    float dt = 1e-2f;
    {
      auto physics_task = profiler.StartScopedTask("[Physics] PreStep", legit::Colors::sunFlower);
      PreStep(
        almost::ParticleGroup::Get(reg),
        almost::LinkGroup::Get(reg),
        almost::TriangleGroup::Get(reg));
    }
    static bool use_simple_solver = true;
    ImGui::Checkbox("Use simple solver", &use_simple_solver);

    auto particle_group = almost::ParticleGroup::Get(reg);
    auto link_group = almost::LinkGroup::Get(reg);
    
    if(use_simple_solver)
    {
      {
        auto physics_task = profiler.StartScopedTask("[Physics] Links", legit::Colors::orange);
        for(int i = 0; i < 10; i++)
        {
          ProjectLinks(
            particle_group.raw<ParticleComponent>(), particle_group.raw<MassComponent>(), particle_group.size(),
            link_group.raw<LinkComponent>(), link_group.raw<LinkIndexComponent>(), link_group.size());
        }
      }
      {
        auto physics_task = profiler.StartScopedTask("[Physics] Integration", legit::Colors::orange);

        auto particle_group = almost::ParticleGroup::Get(reg);
        IntegrateParticles(particle_group.raw<ParticleComponent>(), particle_group.size(), dt);
      }
    }else
    {
      auto particle_group = almost::ParticleGroup::Get(reg);
      auto constraint_graph = ConstraintGraph(
         link_group.raw<LinkIndexComponent>(), link_group.size(),
         particle_group.size()
      );
      std::vector<glm::vec2> inertial_positions;
      inertial_positions.resize(particle_group.size());
      GetParticleInertialPositions(particle_group.raw<ParticleComponent>(), inertial_positions.data(), particle_group.size(), dt);
      AssignInertialPositions(particle_group.raw<ParticleComponent>(), inertial_positions.data(), particle_group.size());
      for(int i = 0; i < 10; i++)
      {
        ProjectVBD(
          particle_group.raw<ParticleComponent>(),
          particle_group.raw<MassComponent>(),
          inertial_positions.data(),
          particle_group.size(),
          link_group.raw<LinkIndexComponent>(),
          link_group.raw<LinkComponent>(),
          link_group.size(),
          constraint_graph,
          dt);
      }
      //AssignInertialPositions(particle_group.raw<ParticleComponent>(), inertial_positions.data(), particle_group.size());
      int p = 1;
    }
  }
}