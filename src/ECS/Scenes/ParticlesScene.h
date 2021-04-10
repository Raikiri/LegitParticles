
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

void CreateGround(entt::registry& reg)
{
  reg.group<almost::ParticleComponent, almost::ParticleIndexComponent>();

  entt::entity groundEntity = CreateParticle(reg, glm::vec2(0.0f, 0.0f), true, false);
}
void CreateClothPhysicsMesh(entt::registry& reg)
{
  reg.group<almost::ParticleComponent, almost::ParticleIndexComponent>();
  glm::vec2 minPoint = glm::vec2(-200.0f, -200.0f);
  glm::vec2 maxPoint = glm::vec2(200.0f, 200.0f);

  glm::ivec2 count = { 40, 40 };
  glm::ivec2 index;
  std::vector<entt::entity> particleEntities;
  auto c = reg.view<almost::ParticleComponent, almost::ParticleIndexComponent>();
  particleEntities.resize(count.x * count.y);
  for (index.y = 0; index.y < count.y; index.y++)
  {
    for (index.x = 0; index.x < count.x; index.x++)
    {
      bool isTopRow = index.y == count.y - 1;
      particleEntities[index.x + index.y * count.x] = CreateParticle(reg, glm::mix(minPoint, maxPoint, glm::vec2(index) / (glm::vec2(count) - glm::vec2(0)))/* + glm::vec2(rand() % 5, rand() % 5)*/, isTopRow, isTopRow);
    }
  }

  reg.group<almost::LinkComponent, almost::LinkIndexComponent>();
  std::vector<entt::entity> linkEntities;
  //particleLinks.resize((count.x - 1) * (count.y - 1));
  for (index.y = 0; index.y < count.y; index.y++)
  {
    for (index.x = 0; index.x < count.x; index.x++)
    {
      entt::entity currEntity = particleEntities[index.x + index.y * count.x];
      if (index.x + 1 < count.x && index.y + 1 < count.y)
      {
        entt::entity nextEntity = particleEntities[index.x + 1 + index.y * count.x];

        linkEntities.push_back(CreateLink(reg, currEntity, nextEntity));
      }
      if (index.y + 1 < count.y)
      {
        entt::entity nextEntity = particleEntities[index.x + (index.y + 1) * count.x];

        linkEntities.push_back(CreateLink(reg, currEntity, nextEntity));
      }
    }
  }
}