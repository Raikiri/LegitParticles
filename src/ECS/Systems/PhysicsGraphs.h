#pragma once
#include "../Components/LinkComponent.h"
#include <vector>

namespace almost
{
  struct ParticleInfo
  {
    size_t links_start = 0;
    size_t links_count = 0;
  };
  
  struct ConstraintGraph
  {
    ConstraintGraph(LinkIndexComponent *link_indices, size_t links_count, size_t particles_count);
    
    using LinkIdx = size_t;
    using ParticleIdx = size_t;
    std::vector<ParticleInfo> particle_infos;
    std::vector<LinkIdx> particle_links;
  };
  
}