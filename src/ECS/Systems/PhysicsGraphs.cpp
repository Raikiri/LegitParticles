#include "PhysicsGraphs.h"
#include <assert.h>

namespace almost
{
  ConstraintGraph::ConstraintGraph(LinkIndexComponent *link_indices, size_t links_count, size_t particles_count)
  {
    particle_links.resize(links_count * 2);
    particle_infos.resize(particles_count);

    for(size_t phase = 0; phase < 2; phase++)
    {
      for(size_t link_idx = 0; link_idx < links_count; link_idx++)
      {
        auto link = link_indices[link_idx];
        for(size_t particle_num = 0; particle_num < 2; particle_num++)
        {
          ParticleIdx particle_idx = link.indices[particle_num];
          assert(particle_idx < particles_count);
          if(phase == 1)
          {
            particle_links[particle_infos[particle_idx].links_start + particle_infos[particle_idx].links_count] = link_idx;
          }
          particle_infos[particle_idx].links_count++;
        }
      }
      if(phase == 0)
      {
        size_t curr_links_offset = 0;
        for(size_t particle_idx = 0; particle_idx < particles_count; particle_idx++)
        {
          particle_infos[particle_idx].links_start = curr_links_offset;
          curr_links_offset += particle_infos[particle_idx].links_count;
          particle_infos[particle_idx].links_count = 0;
        }
        assert(curr_links_offset == links_count * 2);
      }
    }
  }
}