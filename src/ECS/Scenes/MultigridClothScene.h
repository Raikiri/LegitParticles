

entt::entity CreateParticle(entt::registry& reg, glm::vec2 pos, bool isFixed, bool isDraggable)
{
  entt::entity particleEntity = reg.create();
  almost::ParticleComponent particleComponent;
  particleComponent.pos = pos;
#if defined(POSITION_BASED)
  particleComponent.prevPos = pos;
#else
  particleComponent.velocity = glm::vec2(0.0f, 0.0f);
#endif
  particleComponent.acceleration = glm::vec2(0.0f);
  reg.emplace<almost::ParticleComponent>(particleEntity, particleComponent);

  almost::ParticleIndexComponent particleIndexComponent;
  particleIndexComponent.index = -1;
  reg.emplace<almost::ParticleIndexComponent>(particleEntity, particleIndexComponent);


  almost::MassComponent massComponent;
  massComponent.usesGravity = !isFixed;
  massComponent.invMass = isFixed ? 0.0f : 1000.0f;
  reg.emplace<almost::MassComponent>(particleEntity, massComponent);

  almost::DefPosComponent defPosComponent;
  defPosComponent.defPos = pos;
  defPosComponent.isDraggable = isDraggable;
  reg.emplace<almost::DefPosComponent>(particleEntity, defPosComponent);

  reg.emplace<almost::CoarseMultigridComponent>(particleEntity);
  reg.emplace<almost::FineMultigridComponent>(particleEntity);


  return particleEntity;
}

entt::entity CreateLink(entt::registry& reg, entt::entity particle0, entt::entity particle1)
{
  entt::entity linkEntity = reg.create();

  almost::LinkComponent linkComponent;
  linkComponent.entities[0] = particle0;
  linkComponent.entities[1] = particle1;
  linkComponent.defLength = glm::length(reg.get<almost::ParticleComponent>(particle0).pos - reg.get<almost::ParticleComponent>(particle1).pos);
  reg.emplace<almost::LinkComponent>(linkEntity, linkComponent);

  almost::LinkIndexComponent linkIndexComponent;
  linkIndexComponent.indices[0] = -1;
  linkIndexComponent.indices[1] = -1;
  reg.emplace<almost::LinkIndexComponent>(linkEntity, linkIndexComponent);

  return linkEntity;
}

struct TriangleTensorData
{
  almost::StaticTensor<glm::vec2, almost::StrainDynamics::ParticleDim> worldPositions;
  almost::StaticTensor<float, almost::StrainDynamics::ParticleDim> invMasses;
};
TriangleTensorData MakeTriangleTensorData(entt::registry& reg, std::array<entt::entity, 3> particleEntities)
{
  TriangleTensorData data;
  for (auto m : almost::StrainDynamics::ParticleDim())
  {
    data.worldPositions.Get(m) = reg.get<almost::ParticleComponent>(particleEntities[m.value]).pos;
    data.invMasses.Get(m) = reg.get<almost::MassComponent>(particleEntities[m.value]).invMass;
  }
  return data;
}

entt::entity CreateTriangle(entt::registry& reg, entt::entity particle0, entt::entity particle1, entt::entity particle2)
{
  entt::entity triangleEntity = reg.create();

  almost::TriangleComponent triangleComponent;
  triangleComponent.entities[0] = particle0;
  triangleComponent.entities[1] = particle1;
  triangleComponent.entities[2] = particle2;
  TriangleTensorData data = MakeTriangleTensorData(reg, { particle0 , particle1, particle2});
  triangleComponent.uvFromRef = almost::StrainDynamics::BuildUvFromRefMatrix(data.worldPositions);
  reg.emplace<almost::TriangleComponent>(triangleEntity, triangleComponent);

  almost::TriangleIndexComponent triangleIndexComponent;
  triangleIndexComponent.indices[0] = -1;
  triangleIndexComponent.indices[1] = -1;
  triangleIndexComponent.indices[2] = -1;
  reg.emplace<almost::TriangleIndexComponent>(triangleEntity, triangleIndexComponent);


  return triangleEntity;
}

void CreateGround(entt::registry& reg)
{
  reg.group<almost::ParticleComponent, almost::ParticleIndexComponent>();

  entt::entity groundEntity = CreateParticle(reg, glm::vec2(0.0f, 0.0f), true, false);
}
struct MultigridLayer
{
  glm::vec2 minPoint, maxPoint;
  glm::ivec2 count;

  std::vector<entt::entity> linkEntities;
  std::vector<entt::entity> triangleEntities;
  std::vector<entt::entity> particleEntities;
};

MultigridLayer CreateMultigridLayer(entt::registry &reg, glm::vec2 minPoint, glm::vec2 maxPoint, glm::ivec2 count)
{

  MultigridLayer layer;
  layer.minPoint = minPoint;
  layer.maxPoint = maxPoint;
  layer.count = count;

  /*glm::vec2 minPoint = glm::vec2(-200.0f, -200.0f);
  glm::vec2 maxPoint = glm::vec2(200.0f, 200.0f);

  glm::ivec2 count = { 40, 40 };*/
  glm::ivec2 index;
  layer.particleEntities.resize(count.x * count.y);
  for (index.y = 0; index.y < count.y; index.y++)
  {
    for (index.x = 0; index.x < count.x; index.x++)
    {
      bool isTopRow = index.y == count.y - 1;
      bool isFixed = (index.y == count.y - 1) && (index.x == 0 || index.x == count.x - 1);
      layer.particleEntities[index.x + index.y * count.x] = CreateParticle(reg, glm::mix(minPoint, maxPoint, glm::vec2(index) / (glm::vec2(count) - glm::vec2(1)))/* + glm::vec2(rand() % 5, rand() % 5)*/, isFixed, isFixed);
    }
  }

  reg.group<almost::LinkComponent, almost::LinkIndexComponent >();
  reg.group<almost::TriangleComponent, almost::TriangleIndexComponent >();

  //particleLinks.resize((count.x - 1) * (count.y - 1));
  for (index.y = 0; index.y < count.y; index.y++)
  {
    for (index.x = 0; index.x < count.x; index.x++)
    {
      entt::entity currEntity = layer.particleEntities[index.x + index.y * count.x];
      if (index.x + 1 < count.x/* && index.y + 1 < count.y*/)
      {
        entt::entity nextEntity = layer.particleEntities[index.x + 1 + index.y * count.x];

        layer.linkEntities.push_back(CreateLink(reg, currEntity, nextEntity));
      }

      if (index.y + 1 < count.y)
      {
        entt::entity nextEntity = layer.particleEntities[index.x + (index.y + 1) * count.x];

        layer.linkEntities.push_back(CreateLink(reg, currEntity, nextEntity));
      }

      if (index.x + 1 < count.x && index.y + 1 < count.y)
      {
        entt::entity e00 = layer.particleEntities[index.x + index.y * count.x];
        entt::entity e10 = layer.particleEntities[index.x + 1 + index.y * count.x];
        entt::entity e11 = layer.particleEntities[index.x + 1 + (index.y + 1) * count.x];
        entt::entity e01 = layer.particleEntities[index.x + (index.y + 1) * count.x];

        layer.triangleEntities.push_back(CreateTriangle(reg, e00, e10, e11));
        layer.triangleEntities.push_back(CreateTriangle(reg, e00, e11, e01));
      }
    }
  }
  return layer;
}

struct InterpPoint
{
  glm::ivec2 index;
  float weight;
};
std::array<InterpPoint, 4> GetInterpPoints4(glm::ivec2 minIndex, glm::ivec2 maxIndex, glm::vec2 ratio)
{
  std::array<InterpPoint, 4> interpPoints;
  interpPoints[0].index = { minIndex.x, minIndex.y };
  interpPoints[0].weight = (1.0f - ratio.x) * (1.0f - ratio.y);

  interpPoints[1].index = { maxIndex.x, minIndex.y };
  interpPoints[1].weight = ratio.x * (1.0f - ratio.y);

  interpPoints[2].index = { maxIndex.x, maxIndex.y };
  interpPoints[2].weight = ratio.x * ratio.y;

  interpPoints[3].index = { minIndex.x, maxIndex.y };
  interpPoints[3].weight = (1.0f - ratio.x) * ratio.y;
  return interpPoints;
}

entt::entity GetParticleEntity(const MultigridLayer &layer, glm::ivec2 index)
{
  return layer.particleEntities[index.x + index.y * layer.count.x];
}

void AddFineInfluences(almost::ParticleGroup::Type fineParticles, entt::entity fineParticleEntity, almost::ParticleGroup::Type coarseParticles, entt::entity coarseParticleEntity, float weight)
{
  {
    almost::Influence influence;
    influence.particleEntity = coarseParticleEntity;
    influence.weight = weight;
    fineParticles.get<almost::FineMultigridComponent>(fineParticleEntity).influences.push_back(influence);
  }
  /*{
    almost::Influence influence;
    influence.particleEntity = fineParticleEntity;
    influence.weight = weight;
    coarseParticles.get<almost::CoarseMultigridComponent>(coarseParticleEntity).influences.push_back(influence);
  }*/
}


void AddCoarseInfluences(almost::ParticleGroup::Type fineParticles, entt::entity fineParticleEntity, almost::ParticleGroup::Type coarseParticles, entt::entity coarseParticleEntity, float weight)
{
  /*{
    almost::Influence influence;
    influence.particleEntity = coarseParticleEntity;
    influence.weight = weight;
    fineParticles.get<almost::FineMultigridComponent>(fineParticleEntity).influences.push_back(influence);
  }*/
  {
    almost::Influence influence;
    influence.particleEntity = fineParticleEntity;
    influence.weight = weight;
    coarseParticles.get<almost::CoarseMultigridComponent>(coarseParticleEntity).influences.push_back(influence);
  }
}

void NormalizeWeights(std::vector<almost::Influence> &influences)
{
  float totalWeight = 1e-7f;
  for (auto& influence : influences)
  {
    totalWeight += influence.weight;
  }
  for (auto& influence : influences)
  {
    influence.weight /= totalWeight;
  }
}
void StitchMultigridMeshes(const MultigridLayer& fineLayer, almost::ParticleGroup::Type fineParticles, const MultigridLayer& coarseLayer, almost::ParticleGroup::Type coarseParticles)
{
  glm::vec2 fineStep = glm::vec2(1.0f) / (glm::vec2(fineLayer.count) - glm::vec2(1.0f));;
  glm::vec2 coarseStep = glm::vec2(1.0f) / (glm::vec2(coarseLayer.count) - glm::vec2(1.0f));
  glm::ivec2 fineNodeIndex;
  for (fineNodeIndex.y = 0; fineNodeIndex.y < fineLayer.count.y; fineNodeIndex.y++)
  {
    for (fineNodeIndex.x = 0; fineNodeIndex.x < fineLayer.count.x; fineNodeIndex.x++)
    {
      glm::vec2 localCoord = fineStep * glm::vec2(fineNodeIndex);
      glm::vec2 coarseCellIndex2f = localCoord / coarseStep;
      glm::vec2 ratio = glm::fract(coarseCellIndex2f);
      glm::ivec2 minIndex = glm::ivec2(glm::floor(coarseCellIndex2f));
      glm::ivec2 maxIndex = minIndex + glm::ivec2(1);

      minIndex = glm::clamp(minIndex, glm::ivec2(0), coarseLayer.count - glm::ivec2(1));
      maxIndex = glm::clamp(maxIndex, glm::ivec2(0), coarseLayer.count - glm::ivec2(1));
      auto interpPoints = GetInterpPoints4(minIndex, maxIndex, ratio);
      for (auto interpPoint : interpPoints)
      {
        entt::entity fineParticleEntity = GetParticleEntity(fineLayer, fineNodeIndex);
        entt::entity coarseParticleEntity = GetParticleEntity(coarseLayer, interpPoint.index);
        if(abs(interpPoint.weight) > 1e-5f)
        AddFineInfluences(fineParticles, fineParticleEntity, coarseParticles, coarseParticleEntity, interpPoint.weight);
      }
    }
  }

  glm::ivec2 mult = (fineLayer.count - glm::ivec2(1)) / (coarseLayer.count - glm::ivec2(1));

  glm::ivec2 coarseNodeIndex;
  for (coarseNodeIndex.y = 0; coarseNodeIndex.y < coarseLayer.count.y; coarseNodeIndex.y++)
  {
    for (coarseNodeIndex.x = 0; coarseNodeIndex.x < coarseLayer.count.x; coarseNodeIndex.x++)
    {
      glm::ivec2 closestFineNodeIndex = coarseNodeIndex * mult;

      entt::entity fineParticleEntity = GetParticleEntity(fineLayer, closestFineNodeIndex);
      entt::entity coarseParticleEntity = GetParticleEntity(coarseLayer, coarseNodeIndex);
      AddCoarseInfluences(fineParticles, fineParticleEntity, coarseParticles, coarseParticleEntity, 1.0f);
    }
  }

  auto fineMultigridComponents = fineParticles.raw<almost::FineMultigridComponent>();
  for (size_t particleIndex = 0; particleIndex < fineParticles.size(); particleIndex++)
  {
    NormalizeWeights(fineMultigridComponents[particleIndex].influences);
  }

  auto coarseMultigridComponents = coarseParticles.raw<almost::CoarseMultigridComponent>();
  for (size_t particleIndex = 0; particleIndex < coarseParticles.size(); particleIndex++)
  {
    NormalizeWeights(coarseMultigridComponents[particleIndex].influences);
  }
}

void CreateMultigridPhysicsMesh(std::vector<entt::registry> &regLayers, glm::vec2 minPoint, glm::vec2 maxPoint, glm::ivec2 baseSize)
{
  int mult = 2;
  glm::ivec2 currSize = baseSize;
  MultigridLayer prevLayer;
  for (size_t layerIndex = 0; layerIndex < regLayers.size(); layerIndex++)
  {
    auto currLayer = CreateMultigridLayer(regLayers[layerIndex], minPoint, maxPoint, currSize);
    if (layerIndex > 0)
    {
      StitchMultigridMeshes(prevLayer, almost::ParticleGroup::Get(regLayers[layerIndex - 1]), currLayer, almost::ParticleGroup::Get(regLayers[layerIndex]));
    }
    prevLayer = currLayer;
    currSize = (currSize - glm::ivec2(1)) / mult + glm::ivec2(1);
  }
}