# LegitParticles
![multiplatform build](https://github.com/Raikiri/LegitParticles/actions/workflows/cmake-multi-platform.yml/badge.svg)

![image](https://github.com/user-attachments/assets/34713cb3-e3cb-454e-b40f-ea8d4ea15a85)

This is an AVBD implementation faithful to the original paper by Chris Giles, Elie Diaz and Cem Yuskel (link in the description). The demo only link constraints for particles, it's only marginally optimized and runs CPU-side. Both PBD and AVBD solvers use 10 iterations and run in comparable time.

# Impression on AVBD


If you're not familiar with simulating constraint physics, it's particularly challenging to get stiff enough materials such as fabric and ropes that don't stretch. In the demo I show that AVBD without warmstarting is slightly stiffer than PBD. With warmstarting, AVBD is substantially stiffer than PBD.


Unfortunately, the algorithm has quite a few magical parameters that I had to carefully tweak, and they have pretty nasty failure modes such as unstable explosions and weird rigid-like behaviour. Some of the most sensitive parameters are scene dependent not dimensionless:
```cpp
  float beta = 1.0f;
  float min_stiffness = 1e1f;
  float max_stiffness = 1e5f;
  float max_lambda = 1e4f;
```
The paper describes them as if they work in a wide range of values, but in my experience they're quite fiddly and have nasty failure modes.

# Code sample
Probably the most important part of the code (the one doing the actual optimization/integration) is:
```cpp
void ProjectVBD(
    ParticleComponent* particle_components,
    const MassComponent* mass_components,
    const glm::vec2* inertial_positions,
    size_t particles_count,
    const LinkIndexComponent *link_indices,
    const LinkComponent *links,
    size_t links_count,
    const ConstraintGraph &constraint_graph,
    float dt)
  {
    for(size_t particle_idx = 0; particle_idx < particles_count; particle_idx++)
    {
      auto &particle = particle_components[particle_idx];
      auto inv_mass = mass_components[particle_idx].invMass;
      auto inertial_pos = inertial_positions[particle_idx];

      if(inv_mass > 0.0f)
      {
        EnergyDerivatives constraint_derivatives = GetParticleConstraintDerivatives(
          particle.pos,
          particle_idx,
          particle_components,
          link_indices,
          links,
          constraint_graph
        );

        float mdt2 = 1.0f / std::max(1e-7f, dt * dt * inv_mass);
        glm::vec2 f = -constraint_derivatives.grad;
        glm::vec2 total_forces = -(particle.pos - inertial_pos) * mdt2 + f;
        glm::mat2 total_hessian = glm::mat2(mdt2) + constraint_derivatives.hessian;

        if(glm::determinant(total_hessian) > 0.0f)
        {
          glm::vec2 delta = glm::inverse(total_hessian) * total_forces;
          particle.pos += delta;
        }
      }else
      {
        particle.pos = inertial_pos;
      }
    }
  }
```

# Libraries used
LegitVulkan (rendering)

LegitProfiler (profiler UI)

entt (ECS implementation with a memory layout convenient for physics)

imgui (UI)

glfw (window management)

glm (vector maths)

spirv-cross (reflection)
