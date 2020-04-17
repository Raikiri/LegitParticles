#include "../Context/PhysicsData.h"
#include "../Context/MeshRendererData.h"
#include "../Components/ParticleComponent.h"
#include "../Components/LinkComponent.h"
#include <entity\registry.hpp>

#include "../../Utils/GroupArg.h"
#include "../../Maths/AlgebraicMultigridSolver.h"
namespace almost
{
  const glm::vec4 favBlue = glm::vec4(glm::pow(glm::vec3(0.0, 0.5, 0.75), glm::vec3(2.2f)), 1.0f);
  const glm::vec4 favOrange = glm::vec4(glm::pow(glm::vec3(0.8f, 0.4f, 0.1f), glm::vec3(2.2f)), 1.0f);


  PhysicsData InitPhysicsData()
  {
    PhysicsData physicsData;
    return physicsData;
  }

  glm::vec4 Color32(glm::uint32_t col)
  {
    return glm::vec4(
      ((col & 0x000000ff) >> (0 * 8)) / 255.0f,
      ((col & 0x0000ff00) >> (1 * 8)) / 255.0f,
      ((col & 0x00ff0000) >> (2 * 8)) / 255.0f,
      ((col & 0xff000000) >> (3 * 8)) / 255.0f);
  }

  glm::vec4 LinearColor32(glm::uint32_t col)
  {
    return glm::pow(Color32(col), glm::vec4(2.2f));
  }


  void SubmitLine(almost::MeshRendererData& meshRendererData, glm::vec2 point0, glm::vec2 point1, float width, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    glm::vec2 dir = glm::normalize(point1 - point0);
    glm::vec2 tangent = glm::vec2(-dir.y, dir.x);

    meshRendererData.vertexData[vertexOffset + 0] = { glm::vec3(point0 + tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 1] = { glm::vec3(point0 - tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 2] = { glm::vec3(point1 - tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.vertexData[vertexOffset + 3] = { glm::vec3(point1 + tangent * width * 0.5f, 0.0f), glm::vec2(0.0f, 0.0f), color };
    meshRendererData.indexData[indexOffset + 0] = glm::uint32_t(vertexOffset + 0);
    meshRendererData.indexData[indexOffset + 1] = glm::uint32_t(vertexOffset + 1);
    meshRendererData.indexData[indexOffset + 2] = glm::uint32_t(vertexOffset + 2);

    meshRendererData.indexData[indexOffset + 3] = glm::uint32_t(vertexOffset + 0);
    meshRendererData.indexData[indexOffset + 4] = glm::uint32_t(vertexOffset + 2);
    meshRendererData.indexData[indexOffset + 5] = glm::uint32_t(vertexOffset + 3);
    meshRendererData.verticesCount += 4;
    meshRendererData.indicesCount += 6;
  }

  void SubmitCircle(almost::MeshRendererData& meshRendererData, glm::vec2 center, float radius, size_t sectorsCount, glm::vec4 color)
  {
    size_t vertexOffset = meshRendererData.verticesCount;
    size_t indexOffset = meshRendererData.indicesCount;

    for(size_t sectorIndex = 0; sectorIndex < sectorsCount; sectorIndex++)
    {
      float ang = (float(sectorIndex) / float(sectorsCount)) * 2.0f * glm::pi<float>();
      meshRendererData.vertexData[vertexOffset + sectorIndex] = { glm::vec3(center + glm::vec2(cos(ang), sin(ang)) * radius, 0.0f), glm::vec2(0.0f, 0.0f), color };

      if(sectorIndex < sectorsCount - 2)
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 0] = glm::uint32_t(vertexOffset);
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 1] = glm::uint32_t(vertexOffset + sectorIndex + 1);
      meshRendererData.indexData[indexOffset + sectorIndex * 3 + 2] = glm::uint32_t(vertexOffset + sectorIndex + 2);
    }
    meshRendererData.verticesCount += sectorsCount;
    meshRendererData.indicesCount += (sectorsCount - 2) * 3;

  }

  using ParticleGroup = entt::group_type<almost::ParticleComponent, almost::ParticleIndexComponent, almost::MassComponent, almost::DefPosComponent>;
  using LinkGroup = entt::group_type<almost::LinkComponent, almost::LinkIndexComponent>;

  void ProcessPhysicsControls(
    WindowData &windowData,
    ParticleGroup particles,
    InputData &inputData,
    CameraData &cameraData)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();
    DefPosComponent* defPosComponents = particles.raw<DefPosComponent>();

    if (glfwGetKey(windowData.window, GLFW_KEY_SPACE))
    {
      for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
      {
        DefPosComponent& defPosComponent = defPosComponents[particleIndex];
        ParticleComponent& particleComponent = particleComponents[particleIndex];
        particleComponent.pos = defPosComponent.defPos;
        #if defined(POSITION_BASED)
          particleComponent.prevPos = particleComponent.pos;
        #else
          particleComponent.velocity = glm::vec2(0.0f);
        #endif
      }
    }

    glm::vec2 camRightVec = glm::vec2(cos(cameraData.ang), sin(cameraData.ang));
    glm::vec2 camUpVec = glm::vec2(-camRightVec.y, camRightVec.x);
    bool isControlled = glfwGetMouseButton(windowData.window, GLFW_MOUSE_BUTTON_1);
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      DefPosComponent& defPosComponent = defPosComponents[particleIndex];
      MassComponent& massComponent = massComponents[particleIndex];
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      if (!massComponent.usesGravity)
      {
        glm::vec2 dstPos = inputData.worldMousePos + camRightVec * defPosComponent.defPos.x + camUpVec * defPosComponent.defPos.y * 0.0f;
        #if defined(POSITION_BASED)
          particleComponent.pos = isControlled ? dstPos : particleComponent.pos;
        #else
          particleComponent.velocity = isControlled ? ((dstPos - particleComponent.pos) * 100.0f) : glm::vec2(0.0f);
        #endif
      }
    }
  }

  void SolveLinksNonlinearGauss(ParticleGroup particles, LinkGroup links)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    LinkComponent* linkComponents = links.raw<LinkComponent>();

    for (int i = 0; i < 10; i++)
    {
      for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
      {
        size_t particleIndex0 = linkIndexComponents[linkIndex].indices[0];
        MassComponent& massComponent0 = massComponents[particleIndex0];
        ParticleComponent& particleComponent0 = particleComponents[particleIndex0];

        size_t particleIndex1 = linkIndexComponents[linkIndex].indices[1];
        MassComponent& massComponent1 = massComponents[particleIndex1];
        ParticleComponent& particleComponent1 = particleComponents[particleIndex1];

        glm::vec2 delta = glm::normalize(particleComponent0.pos - particleComponent1.pos);
        float compInvMass = 1.0f / (massComponent0.invMass + massComponent1.invMass + 1e-7f);

        float deltaAcceleration = glm::dot(delta, particleComponent0.acceleration - particleComponent1.acceleration);
        if (fabs(deltaAcceleration) > 1e-2f)
        {
          int pp = 1;
        }
        //float deltaVelocity = glm::dot(delta, particleComponent0.velocity - particleComponent1.velocity);
        float deltaPos = glm::dot(delta, particleComponent0.pos - particleComponent1.pos) - linkComponents[linkIndex].defLength;

        float lambdaAcceleration = deltaAcceleration * compInvMass;
        //float lambdaVelocity = deltaVelocity * compInvMass;
        float lambdaPos = deltaPos * compInvMass;

        particleComponent0.pos += -delta * lambdaPos * massComponent0.invMass;
        particleComponent1.pos +=  delta * lambdaPos * massComponent1.invMass;
        /*particleComponent0.velocity += -delta * lambdaVelocity * massComponent0.invMass;
        particleComponent1.velocity += delta * lambdaVelocity * massComponent1.invMass;
        particleComponent0.acceleration += -delta * lambdaAcceleration * massComponent0.invMass;
        particleComponent1.acceleration += delta * lambdaAcceleration * massComponent1.invMass;*/
      }
    }
  }



  void BuildSparseJacobian(ParticleGroup particles, LinkGroup links,
    PhysicsData::JacobianMatrix& jacobianMatrix, float* positionRightSideVector, float* velocityRightSideVector, float *accelerationRightSide)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    LinkComponent* linkComponents = links.raw<LinkComponent>();

    for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
    {
      size_t particleIndex0 = linkIndexComponents[linkIndex].indices[0];
      MassComponent& massComponent0 = massComponents[particleIndex0];
      ParticleComponent& particleComponent0 = particleComponents[particleIndex0];

      size_t particleIndex1 = linkIndexComponents[linkIndex].indices[1];
      MassComponent& massComponent1 = massComponents[particleIndex1];
      ParticleComponent& particleComponent1 = particleComponents[particleIndex1];

      glm::vec2 delta = glm::normalize(particleComponent0.pos - particleComponent1.pos);

      float deltaAcceleration = -glm::dot(delta, particleComponent0.acceleration - particleComponent1.acceleration);
      #if defined(POSITION_BASED)
        float deltaVelocity = 0.0f;
      #else
        float deltaVelocity = -glm::dot(delta, particleComponent0.velocity - particleComponent1.velocity);
      #endif
      float deltaPos = -(glm::dot(delta, particleComponent0.pos - particleComponent1.pos) - linkComponents[linkIndex].defLength);

      {
        PhysicsData::ParticleIndex columnIndices[2];
        columnIndices[0] = particleIndex0;
        columnIndices[1] = particleIndex1;

        Vector2f jacobianValues[2];
        jacobianValues[0] = {  delta };
        jacobianValues[1] = { -delta };
        jacobianMatrix.AddRow(linkIndex, columnIndices, jacobianValues, 2);
      }
      {
        positionRightSideVector[linkIndex] = deltaPos;
      }
      {
        velocityRightSideVector[linkIndex] = deltaVelocity;
      }
      {
        accelerationRightSide[linkIndex] = deltaAcceleration;
      }
    }
    //jacobianMatrix.Finalize();
    assert(jacobianMatrix.CheckSortedIndices());
  }

  void BuildSparseMassMatrix(ParticleGroup particles, PhysicsData::MassMatrix& massMatrix)
  {
    size_t particlesCount = particles.size();
    MassComponent* massComponents = particles.raw<MassComponent>();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      massMatrix.AddRow(particleIndex, &particleIndex, &(massComponents[particleIndex].invMass), 1);
    }
    //massMatrix.Finalize();
    assert(massMatrix.CheckSortedIndices());
  }

  void BuildSparseVectorFromDense(float* values, size_t count, PhysicsData::RightSideMatrix& sparseVector)
  {
    for (size_t rowIndex = 0; rowIndex < count; rowIndex++)
    {
      char columnIndex = 0;
      float value = values[rowIndex];
      sparseVector.AddRow(rowIndex, &columnIndex, &value, 1);
    }
  }

  void GaussSolveSparseSystem(const PhysicsData::SystemMatrix& systemMatrix, const float* rightSide, float* lambdas, size_t iterationsCount)
  {
    std::vector<PhysicsData::LinkIndex> adjacentLinkIndices;
    std::vector<float> adjacentLinkWeights;
    size_t linksCount = systemMatrix.GetRowsCount();
    for (size_t iterationIndex = 0; iterationIndex < iterationsCount; iterationIndex++)
    {
      for (size_t linkIndex = 0; linkIndex < linksCount; linkIndex++)
      {
        size_t adjacentLinksCount = systemMatrix.GetRowTermsCount(linkIndex);
        adjacentLinkIndices.resize(adjacentLinksCount);
        adjacentLinkWeights.resize(adjacentLinksCount);
        systemMatrix.GetRowTerms(linkIndex, adjacentLinkIndices.data(), adjacentLinkWeights.data());
        float remainder = rightSide[linkIndex];
        float ownCoeff = -1.0f;

        for (size_t adjacentLinkNumber = 0; adjacentLinkNumber < adjacentLinksCount; adjacentLinkNumber++)
        {
          size_t adjacentLinkIndex = adjacentLinkIndices[adjacentLinkNumber];
          if (adjacentLinkIndex == linkIndex)
            ownCoeff = adjacentLinkWeights[adjacentLinkNumber];
          else
            remainder -= adjacentLinkWeights[adjacentLinkNumber] * lambdas[adjacentLinkIndex];
        }
        assert(ownCoeff > 0.0f);
        lambdas[linkIndex] = remainder / ownCoeff;
      }
    }
  }
  void ApplyAccelerationDeltas(ParticleGroup particles, const PhysicsData::DeltaMatrix& deltaMatrix)
  {
    size_t rowsCount = deltaMatrix.GetRowsCount();
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(deltaMatrix.GetRowTermsCount(rowNumber) == 1);
      char columnIndex;
      Vector2f delta;
      deltaMatrix.GetRowTerms(rowNumber, &columnIndex, &delta);
      assert(columnIndex == 0);
      particleComponents[rowNumber].acceleration += delta.vec;
    }
  }
  void ApplyVelocityDeltas(ParticleGroup particles, const PhysicsData::DeltaMatrix& deltaMatrix)
  {
    size_t rowsCount = deltaMatrix.GetRowsCount();
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(deltaMatrix.GetRowTermsCount(rowNumber) == 1);
      char columnIndex;
      Vector2f delta;
      deltaMatrix.GetRowTerms(rowNumber, &columnIndex, &delta);
      assert(columnIndex == 0);
      #if defined(POSITION_BASED)
        assert(0);
      #else
        particleComponents[rowNumber].velocity += delta.vec;
      #endif
    }
  }

  void ApplyPositionDeltas(ParticleGroup particles, const PhysicsData::DeltaMatrix& deltaMatrix)
  {
    size_t rowsCount = deltaMatrix.GetRowsCount();
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(deltaMatrix.GetRowTermsCount(rowNumber) == 1);
      char columnIndex;
      Vector2f delta;
      deltaMatrix.GetRowTerms(rowNumber, &columnIndex, &delta);
      assert(columnIndex == 0);
      particleComponents[rowNumber].pos += delta.vec;
    }
  }
  Vector2f operator *(Vector2f left, float right)
  {
    return Vector2f{ left.vec * right };
  }
  Vector2f operator *(float left, Vector2f right)
  {
    return Vector2f{ right.vec * left };
  }
  float operator *(Vector2f left, Vector2f right)
  {
    return glm::dot(left.vec, right.vec);
  }
  Vector2f& operator += (Vector2f &left, Vector2f right)
  {
    left.vec += right.vec;
    return left;
  }
  void SolveLinksMultigrid(ParticleGroup particles, LinkGroup links)
  {
    using ElemType0 = almost::SparseMatrix<size_t, size_t, Vector2f>::Element;
    using ElemType1 = almost::SparseMatrix<size_t, size_t, float>::Element;
    using ElemType2 = almost::SparseMatrix<size_t, char, Vector2f>::Element;
    using StorageType = almost::StackStorage<std::vector<float>, std::vector<size_t>, std::vector<Vector2f>, std::vector<ElemType0>, std::vector<ElemType1>, std::vector<ElemType2>, PhysicsData::SystemMatrix>;
    StorageType stackStorage;


    std::vector<float> accelerationLambdas;
    std::vector<float> velocityLambdas;
    std::vector<float> positionLambdas;

    std::vector<float> positionRightSide;
    std::vector<float> velocityRightSide;
    std::vector<float> accelerationRightSide;

    PhysicsData::JacobianMatrix jacobianMatrix;
    PhysicsData::MassMatrix massMatrix;
    PhysicsData::JacobianTransposedMatrix jacobianTransposed;
    PhysicsData::JacobianTransposedMatrix massJacobianTransposed;
    PhysicsData::SystemMatrix systemMatrix;

    PhysicsData::DeltaMatrix velocityDeltaMatrix;
    PhysicsData::DeltaMatrix positionDeltaMatrix;
    PhysicsData::DeltaMatrix accelerationDeltaMatrix;

    positionRightSide.resize(links.size());
    velocityRightSide.resize(links.size());
    accelerationRightSide.resize(links.size());

    BuildSparseJacobian(particles, links, jacobianMatrix, positionRightSide.data(), velocityRightSide.data(), accelerationRightSide.data());

    BuildSparseMassMatrix(particles, massMatrix);

    jacobianTransposed.BuildFromTransposed(jacobianMatrix, stackStorage);
    massJacobianTransposed.BuildFromSparseProduct(massMatrix, jacobianTransposed, stackStorage);
    systemMatrix.BuildFromSparseProduct(jacobianMatrix, massJacobianTransposed, stackStorage);

    accelerationLambdas.resize(links.size());
    size_t iterationsCount = 5;
    almost::AlgebraicMultigridSolver<PhysicsData::LinkIndex, float> multigridSolver;
    multigridSolver.LoadSystem(systemMatrix, links.size(), stackStorage);

    /*if (0)
    {
      GaussSolveSparseSystem(systemMatrix, accelerationRightSide.data(), accelerationLambdas.data(), iterationsCount);
    }
    else
    {

      multigridSolver.Solve(accelerationRightSide.data(), iterationsCount, stackStorage);

      std::copy_n(multigridSolver.GetValues(), links.size(), accelerationLambdas.data());
    }

    PhysicsData::RightSideMatrix accelerationLambdasSparse;
    BuildSparseVectorFromDense(accelerationLambdas.data(), accelerationLambdas.size(), accelerationLambdasSparse);
    accelerationDeltaMatrix.BuildFromSparseProduct(massJacobianTransposed, accelerationLambdasSparse, stackStorage);
    ApplyAccelerationDeltas(particles, accelerationDeltaMatrix);*/

    #if !defined(POSITION_BASED)
      velocityLambdas.resize(links.size());
      if (0)
      {
        GaussSolveSparseSystem(systemMatrix, velocityRightSide.data(), velocityLambdas.data(), iterationsCount);
      }
      else
      {
        multigridSolver.Solve(velocityRightSide.data(), iterationsCount, stackStorage);

        std::copy_n(multigridSolver.GetValues(), links.size(), velocityLambdas.data());
      }
      PhysicsData::RightSideMatrix velocityLambdasSparse;
      BuildSparseVectorFromDense(velocityLambdas.data(), velocityLambdas.size(), velocityLambdasSparse);
      velocityDeltaMatrix.BuildFromSparseProduct(massJacobianTransposed, velocityLambdasSparse, stackStorage);
      ApplyVelocityDeltas(particles, velocityDeltaMatrix);
    #endif

    positionLambdas.resize(links.size());
    if (0)
    {
      GaussSolveSparseSystem(systemMatrix, positionRightSide.data(), positionLambdas.data(), iterationsCount);
    }
    else
    {
      multigridSolver.Solve(positionRightSide.data(), iterationsCount, stackStorage);

      std::copy_n(multigridSolver.GetValues(), links.size(), positionLambdas.data());
    }

    PhysicsData::RightSideMatrix positionLambdasSparse;
    BuildSparseVectorFromDense(positionLambdas.data(), positionLambdas.size(), positionLambdasSparse);
    positionDeltaMatrix.BuildFromSparseProduct(massJacobianTransposed, positionLambdasSparse, stackStorage);
    ApplyPositionDeltas(particles, positionDeltaMatrix);
  }

  void ProcessPhysics(
    ParticleGroup particles,
    LinkGroup links,
    almost::PhysicsData& physicsData)
  {
    float dt = 1e-2f;
    glm::vec3 gravity = { 0.0f, -10000.0f, 0.0f };
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    ParticleIndexComponent* particleIndicesComponents = particles.raw<ParticleIndexComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      MassComponent& massComponent = massComponents[particleIndex];
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      particleComponent.acceleration = massComponent.usesGravity ? gravity : glm::vec3(0.0f, 0.0f, 0.0f);
      particleIndicesComponents[particleIndex].index = particleIndex;
    }

    LinkComponent* linkComponents = links.raw<LinkComponent>();
    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
    {
      auto particleEntity0 = linkComponents[linkIndex].entities[0];
      auto particleEntity1 = linkComponents[linkIndex].entities[1];
      assert(particles.contains(particleEntity0));
      assert(particles.contains(particleEntity1));
      linkIndexComponents[linkIndex].indices[0] = particles.get<ParticleIndexComponent>(particleEntity0).index;
      linkIndexComponents[linkIndex].indices[1] = particles.get<ParticleIndexComponent>(particleEntity1).index;
    }
    /*for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      particleComponent.velocity += particleComponent.acceleration * dt;
    }*/
    //SolveLinksNonlinearGauss(particles, links);
    SolveLinksMultigrid(particles, links);

    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      ParticleComponent& particleComponent = particleComponents[particleIndex];
      #if defined(POSITION_BASED)
        glm::vec2 tmpPos = particleComponent.pos;
        particleComponent.pos += (particleComponent.pos - particleComponent.prevPos) + particleComponent.acceleration * dt * dt;
        particleComponent.prevPos = tmpPos;
      #else
        particleComponent.velocity += particleComponent.acceleration * dt;
        particleComponent.pos += particleComponent.velocity * dt;
      #endif
    }
  }

  void SubmitParticles(
    ParticleGroup &particles,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    ParticleIndexComponent* particlesIndicesComponents = particles.raw<ParticleIndexComponent>();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      SubmitCircle(meshRendererData, particleComponents[particleIndex].pos, 8.0f, 15, LinearColor32(legit::Colors::belizeHole)/*favBlue*/);
    }
  }
  void SubmitLinks(
    ParticleGroup &particles,
    LinkGroup& links,
    almost::PhysicsData& physicsData,
    almost::MeshRendererData& meshRendererData)
  {
    LinkComponent* linkComponents = links.raw<LinkComponent>();
    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
    {
      auto particleEntity0 = linkComponents[linkIndex].entities[0];
      auto particleEntity1 = linkComponents[linkIndex].entities[1];
      assert(particles.contains(particleEntity0));
      assert(particles.contains(particleEntity1));
      ParticleComponent particleComponent0 = particles.get<ParticleComponent>(particleEntity0);
      ParticleComponent particleComponent1 = particles.get<ParticleComponent>(particleEntity1);
      SubmitLine(meshRendererData, particleComponent0.pos, particleComponent1.pos, 4.0f, LinearColor32(legit::Colors::orange)/*favOrange*/);
    }
  }

  void ProcessPhysicsData(almost::PhysicsData& physicsData, almost::MeshRendererData& meshRendererData)
  {
    /*almost::MeshRendererData::Vertex vertex;
    meshRendererData.vertexData[meshRendererData.verticesCount++] = { glm::vec3(0.0f, 0.00f, 0.0f), glm::vec2(0.0f, 0.0f), LinearColor32(legit::Colors::carrot) };
    meshRendererData.vertexData[meshRendererData.verticesCount++] = { glm::vec3(0.0f, 450.0f, 0.0f), glm::vec2(0.0f, 0.0f), LinearColor32(legit::Colors::carrot) };
    meshRendererData.vertexData[meshRendererData.verticesCount++] = { glm::vec3(450.0f, 0.0f, 0.0f), glm::vec2(0.0f, 0.0f), LinearColor32(legit::Colors::carrot) };
    meshRendererData.indexData[0] = 0;
    meshRendererData.indexData[1] = 1;
    meshRendererData.indexData[2] = 2;
    meshRendererData.indicesCount = 3;*/
    SubmitLine(meshRendererData, glm::vec2(-450.0f, 50.0f), glm::vec2(450.0f, 0.0f), 10.0f, LinearColor32(legit::Colors::carrot));
  }

  void SubmitGrid(almost::PhysicsData& physicsData, almost::MeshRendererData& meshRendererData, almost::CameraData& cameraData, float step, glm::vec4 color)
  {
  }
}