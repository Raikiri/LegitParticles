#pragma once

namespace almost
{
  struct DotProduct
  {
    using ValueType0 = glm::vec2;
    using ValueType1 = glm::vec2;
    using ResultType = float;
    using ScalarType = float;

    static ResultType Mul(glm::vec2 val0, glm::vec2 val1)
    {
      return glm::dot(val0, val1);
    }
  };

  struct VectorScalarProduct
  {
    using ValueType0 = glm::vec2;
    using ValueType1 = float;
    using ResultType = glm::vec2;
    using ScalarType = float;

    static ResultType Mul(glm::vec2 val0, float val1)
    {
      return val0 * val1;
    }
    static ResultType Mul(float val0, glm::vec2 val1)
    {
      return val0 * val1;
    }
  };

  struct TensorVectorProduct
  {
    using ValueType0 = Tensor2f;
    using ValueType1 = glm::vec2;
    using ResultType = glm::vec2;
    using ScalarType = float;

    static ResultType Mul(Tensor2f val0, glm::vec2 val1)
    {
      return val0(val1);
    }
    static glm::vec2 Divide(ResultType val0, Tensor2f val1)
    {
      return val1.GetInverse()(val0);
    }
    static Tensor2f Mul(ScalarType scalar, Tensor2f val)
    {
      return val * scalar;
    }
    static Tensor2f Mul(Tensor2f val, ScalarType scalar)
    {
      return val * scalar;
    }


    static float SpectralRadius(Tensor2f center, Tensor2f val0, Tensor2f val1)
    {
      float rc = sqrt(center.xx * center.xx + center.yy * center.yy);
      float r0 = sqrt(val0.xx * val0.xx + val0.yy * val0.yy);
      float r1 = sqrt(val1.xx * val1.xx + val1.yy * val1.yy);
      return rc / (sqrt(r0) * sqrt(r1));
    }
  };


  using StorageType = almost::StackStorage<
    PhysicsData::SystemMatrix,
    PhysicsData::MassMatrix,
    PhysicsData::JacobianMatrix,
    PhysicsData::JacobianTransposedMatrix,
    PhysicsData::DeltaMatrix,
    PhysicsData::RightSideMatrix,
    std::vector<PhysicsData::DeltaMatrix::Element>,
    std::vector<PhysicsData::SystemMatrix::Element>,
    std::vector<PhysicsData::JacobianTransposedMatrix::Element>,
    std::vector<PhysicsData::ParticleTensorMatrix::Element>,
    std::vector<float>,
    std::vector<size_t>,
    std::vector<glm::vec2> >;

  using ImplicitStorageType = almost::StackStorage<
    PhysicsData::ParticleScalarMatrix,
    std::vector<PhysicsData::ParticleScalarMatrix::Element>,
    std::vector<PhysicsData::ParticleTensorMatrix::Element>,
    std::vector<PhysicsData::DeltaMatrix::Element>,
    std::vector<Tensor2f>,
    std::vector<size_t>,
    std::vector<float>>;


  PhysicsData::JacobianTransposedMatrix massJacobianTransposed;
  PhysicsData::SystemMatrix systemMatrix;
  almost::AlgebraicMultigridSolver<BasicSpace<float>, PhysicsData::JointDimension> multigridSolver;
  StorageType stackStorage;


  almost::AlgebraicMultigridSolver<TensorVectorProduct, PhysicsData::ParticleDimension> implicitMultigridSolver;
  ImplicitStorageType implicitStorage;



  void BuildSparseJacobian(ParticleGroup::Type particles, LinkGroup::Type links,
    PhysicsData::JacobianMatrix& jacobianMatrix, float* positionRightSideVector, float* velocityRightSideVector, float* accelerationRightSide)
  {
    jacobianMatrix.BuildEmpty(PhysicsData::JointDimension{ links.size() }, PhysicsData::ParticleDimension{ particles.size() });
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
        jacobianMatrix.AppendRow(linkIndex);
        jacobianMatrix.AppendTerm(particleIndex0, delta);
        jacobianMatrix.AppendTerm(particleIndex1, -delta);
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

  void BuildSparseMassMatrix(ParticleGroup::Type particles, PhysicsData::MassMatrix& massMatrix)
  {
    massMatrix.BuildEmpty(PhysicsData::ParticleDimension{ particles.size() }, PhysicsData::ParticleDimension{ particles.size() });
    size_t particlesCount = particles.size();
    MassComponent* massComponents = particles.raw<MassComponent>();
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      massMatrix.AppendRow(particleIndex);
      massMatrix.AppendTerm(particleIndex, massComponents[particleIndex].invMass);
    }
    //massMatrix.Finalize();
    assert(massMatrix.CheckSortedIndices());
  }

  template<typename ValueType, typename MatrixType>
  void BuildSparseVectorFromDense(ValueType* values, size_t count, MatrixType& sparseVector)
  {
    sparseVector.BuildEmpty(typename MatrixType::RowDimension{ count }, PhysicsData::OneDimension());
    for (size_t rowIndex = 0; rowIndex < count; rowIndex++)
    {
      char columnIndex = 0;
      float value = values[rowIndex];
      sparseVector.AppendRow(rowIndex);
      sparseVector.AppendTerm(0, value);
    }
  }
  void ApplyAccelerationDeltas(ParticleGroup::Type particles, const PhysicsData::DeltaMatrix& deltaMatrix)
  {
    size_t rowsCount = deltaMatrix.GetRowsCount();
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(deltaMatrix.GetRowTermsCount(rowNumber) == 1);
      auto rowTerms = deltaMatrix.GetRowTerms(rowNumber);
      size_t rowIndex = deltaMatrix.GetRowIndex(rowNumber);
      glm::vec2 delta = rowTerms[0].value;
      assert(rowTerms[0].columnIndex == 0);
      particleComponents[rowIndex].acceleration += delta;
    }
  }
  void ApplyVelocityDeltas(ParticleGroup::Type particles, const PhysicsData::DeltaMatrix& deltaMatrix)
  {
    size_t rowsCount = deltaMatrix.GetRowsCount();
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(deltaMatrix.GetRowTermsCount(rowNumber) == 1);
      auto rowTerms = deltaMatrix.GetRowTerms(rowNumber);
      size_t rowIndex = deltaMatrix.GetRowIndex(rowNumber);
      glm::vec2 delta = rowTerms[0].value;
      assert(rowTerms[0].columnIndex == 0);
#if defined(POSITION_BASED)
      assert(0);
#else
      particleComponents[rowIndex].velocity += delta;
#endif
    }
  }

  void ApplyPositionDeltas(ParticleGroup::Type particles, const PhysicsData::DeltaMatrix& deltaMatrix)
  {
    size_t rowsCount = deltaMatrix.GetRowsCount();
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(deltaMatrix.GetRowTermsCount(rowNumber) == 1);
      auto rowTerms = deltaMatrix.GetRowTerms(rowNumber);
      size_t rowIndex = deltaMatrix.GetRowIndex(rowNumber);
      glm::vec2 delta = rowTerms[0].value;
      assert(rowTerms[0].columnIndex == 0);
      particleComponents[rowIndex].pos += delta;
    }
  }

  void ApplyPositionDenseDeltas(ParticleGroup::Type particles, const glm::vec2* deltaVector)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      glm::vec2 delta = deltaVector[particleIndex];
      particleComponents[particleIndex].pos += delta;
    }
  }
  void ApplyVelocityDenseDeltas(ParticleGroup::Type particles, const glm::vec2* deltaVector)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
#if defined(POSITION_BASED)
      assert(0);
#else
      glm::vec2 delta = deltaVector[particleIndex];
      /*float maxLen = 1000.0f;
      if (glm::length(delta) > maxLen)
        delta = glm::normalize(delta) * maxLen;*/
      particleComponents[particleIndex].velocity += delta;
#endif
    }
  }



  template<typename ValueType, typename MatrixType>
  void BuildDenseVectorFromSparse(ValueType* values, const MatrixType& sparseVector)
  {
    assert(sparseVector.columnDimension.size == 1);

    size_t rowsCount = sparseVector.GetRowsCount();
    std::fill_n(values, sparseVector.rowDimension.size, ValueType(0));
    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(sparseVector.GetRowTermsCount(rowNumber) == 1);
      auto rowTerms = sparseVector.GetRowTerms(rowNumber);
      assert(rowTerms[0].columnIndex == 0);
      values[sparseVector.GetRowIndex(rowNumber)] = rowTerms[0].value;
    }
  }

  void SolveLinksMultigrid(ParticleGroup::Type particles, LinkGroup::Type links, legit::CpuProfiler& profiler)
  {
    auto positionRightSideHandle = stackStorage.GetHandle< std::vector<float> >();
    auto velocityRightSideHandle = stackStorage.GetHandle< std::vector<float> >();
    auto accelerationRightSideHandle = stackStorage.GetHandle< std::vector<float> >();
    auto positionRightSide = positionRightSideHandle.Get();
    auto velocityRightSide = velocityRightSideHandle.Get();
    auto accelerationRightSide = accelerationRightSideHandle.Get();

    positionRightSide.resize(links.size());
    velocityRightSide.resize(links.size());
    accelerationRightSide.resize(links.size());

    size_t iterationsCount = 10;

    {
      auto massMatrixHandle = stackStorage.GetHandle<PhysicsData::MassMatrix>();
      auto& massMatrix = massMatrixHandle.Get();

      auto jacobianMatrixHandle = stackStorage.GetHandle<PhysicsData::JacobianMatrix>();
      auto& jacobianMatrix = jacobianMatrixHandle.Get();

      auto jacobianTransposedMatrixHandle = stackStorage.GetHandle<PhysicsData::JacobianTransposedMatrix>();
      auto& jacobianTransposedMatrix = jacobianTransposedMatrixHandle.Get();

      {
        auto physicsTask = profiler.StartScopedTask("[Physics] Jacobian", legit::Colors::greenSea);
        BuildSparseJacobian(particles, links, jacobianMatrix, positionRightSide.data(), velocityRightSide.data(), accelerationRightSide.data());

      }
      {
        auto physicsTask = profiler.StartScopedTask("[Physics] System", legit::Colors::sunFlower);

        BuildSparseMassMatrix(particles, massMatrix);

        jacobianTransposedMatrix.BuildFromTransposed(jacobianMatrix, stackStorage);
        massJacobianTransposed.BuildFromSparseProduct<VectorScalarProduct>(massMatrix, jacobianTransposedMatrix, stackStorage);
        systemMatrix.BuildFromSparseProduct<DotProduct>(jacobianMatrix, massJacobianTransposed, stackStorage);
        /*massJacobianTransposed.BuildFromDenseProduct<VectorScalarProduct>(massMatrix, jacobianMatrix, stackStorage);
        systemMatrix.BuildFromSparseProduct<DotProduct>(jacobianMatrix, massJacobianTransposed, stackStorage);*/
        int p = 1;
      }
    }

    {
      auto physicsTask = profiler.StartScopedTask("[Physics] Load", legit::Colors::pomegranate);
      multigridSolver.LoadSystem(systemMatrix, stackStorage);
    }

    {
      auto physicsTask = profiler.StartScopedTask("[Physics] Solve", legit::Colors::amethyst);

      /*
      accelerationLambdas.resize(links.size());
      if (0)
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
      auto velocityLambdasSparseHandle = stackStorage.GetHandle<PhysicsData::RightSideMatrix>();
      auto velocityDeltaMatrixHandle = stackStorage.GetHandle<PhysicsData::DeltaMatrix>();

      auto& velocityLambdasSparse = velocityLambdasSparseHandle.Get();
      auto& velocityDeltaMatrix = velocityDeltaMatrixHandle.Get();
      if (0)
      {
        auto velocityLambdasHandle = stackStorage.GetHandle<std::vector<float>>();
        auto& velocityLambdas = velocityLambdasHandle.Get();
        velocityLambdas.resize(links.size());

        GaussSolveSparseSystem(systemMatrix, velocityRightSide.data(), velocityLambdas.data(), iterationsCount);
        BuildSparseVectorFromDense(velocityLambdas.data(), velocityLambdas.size(), velocityLambdasSparse);
      }
      else
      {
        //multigridSolver.Solve(velocityRightSide.data(), iterationsCount, stackStorage);
        //BuildSparseVectorFromDense(multigridSolver.GetValues(), links.size(), velocityLambdasSparse);
      }

      velocityDeltaMatrix.BuildFromSparseProduct<VectorScalarProduct>(massJacobianTransposed, velocityLambdasSparse, stackStorage);
      ApplyVelocityDeltas(particles, velocityDeltaMatrix);
#endif

      {
        auto positionLambdasSparseHandle = stackStorage.GetHandle<PhysicsData::RightSideMatrix>();
        auto positionDeltaMatrixHandle = stackStorage.GetHandle<PhysicsData::DeltaMatrix>();

        auto& positionLambdasSparse = positionLambdasSparseHandle.Get();
        auto& positionDeltaMatrix = positionDeltaMatrixHandle.Get();
        if (1)
        {
          auto positionLambdasHandle = stackStorage.GetHandle<std::vector<float>>();
          auto& positionLambdas = positionLambdasHandle.Get();
          positionLambdas.resize(links.size());

          IterateGaussSeidel<BasicSpace<float>>(systemMatrix, positionRightSide.data(), positionLambdas.data(), 100, implicitStorage);

          BuildSparseVectorFromDense(positionLambdas.data(), positionLambdas.size(), positionLambdasSparse);
        }
        else
        {
          multigridSolver.Solve(positionRightSide.data(), iterationsCount, stackStorage);
          BuildSparseVectorFromDense(multigridSolver.GetValues(), links.size(), positionLambdasSparse);
        }

        positionDeltaMatrix.BuildFromSparseProduct<VectorScalarProduct>(massJacobianTransposed, positionLambdasSparse, stackStorage);
        ApplyPositionDeltas(particles, positionDeltaMatrix);
      }
    }
  }




  struct ConstraintData
  {
    Tensor2f derivative;
    glm::vec2 delta;
  };
  ConstraintData GetLinkDataLinear(glm::vec2 pos0, glm::vec2 pos1, float defLen)
  {

    //C = sqrt((pj - qj)(pj - qj)) - d
    //dC/dpi = (pi - qi)/|p - q|
    //d/drj ri/sqrt = Dij /sqrt - rirj / r^3
    glm::vec2 delta = pos1 - pos0;
    glm::vec2 dir = glm::normalize(delta);

    //linear
    auto unitTensor = Tensor2f(1.0f);
    float len = glm::length(delta);
    ConstraintData linkData;
    float constraint = (len - defLen);
    linkData.delta = dir * constraint;
    linkData.derivative = almost::TensorProduct<float>(dir) * -1.0f + (unitTensor - TensorProduct<float>(dir)) * (-constraint / len);
    return linkData;
  }

  ConstraintData GetLinkDataQuadratic(glm::vec2 pos0, glm::vec2 pos1, float defLen)
  {
    //C=(pos0 - pos1)^2 - d^2
    //dC/dpi = 2Dij(pj - qj)=2(pi - qi)
    //d2C/dpi dqj = 2 Dij
    auto unitTensor = Tensor2f(1.0f);
    glm::vec2 delta = pos1 - pos0;
    ConstraintData linkData;
    float constraint = (glm::dot(delta, delta) - defLen * defLen);
    linkData.delta = delta * 2.0f * constraint;
    linkData.derivative = almost::TensorProduct<float>(delta) * -4.0f + unitTensor * -2.0f * constraint;
    return linkData;
  }
  //https://www.cs.cmu.edu/~baraff/papers/sig98.pdf
  //Large Steps in Cloth Simulation David Baraff Andrew Witkin
  void BuildConstraintMatrices(ParticleGroup::Type particles, LinkGroup::Type links,
    PhysicsData::ParticleTensorMatrix& constraintDerivatives, glm::vec2* constraintDeltas)
  {
    glm::vec2 gravity = glm::vec2(0.0f, -10000.0f);
    constraintDerivatives.BuildEmpty(PhysicsData::ParticleDimension{ particles.size() }, PhysicsData::ParticleDimension{ particles.size() });

    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    LinkIndexComponent* linkIndexComponents = links.raw<LinkIndexComponent>();
    LinkComponent* linkComponents = links.raw<LinkComponent>();

    using Element = PhysicsData::ParticleTensorMatrix::Element;
    std::vector<Element> elements;
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      constraintDeltas[particleIndex] = (massComponents[particleIndex].invMass > 0.0f) ? gravity : glm::vec2(0.0f);
    }
    float linkStiffness = 100.1f;
    for (size_t linkIndex = 0; linkIndex < links.size(); linkIndex++)
    {
      auto& link = linkComponents[linkIndex];
      size_t particleIndices[2];
      glm::vec2 particlePos[2];
      float particleInvMass[2];
      for (auto i = 0; i < 2; i++)
      {
        particleIndices[i] = linkIndexComponents[linkIndex].indices[i];
        particlePos[i] = particleComponents[particleIndices[i]].pos;
        particleInvMass[i] = massComponents[particleIndices[i]].invMass;
      }

      auto linkData = GetLinkDataLinear(particlePos[0], particlePos[1], link.defLength);
      //auto linkData = GetLinkDataQuadratic(particlePos[0], particlePos[1], link.defLength);

      Element e;
      e.rowIndex = particleIndices[0];
      e.columnIndex = particleIndices[1];
      e.value = linkData.derivative * particleInvMass[0] * -linkStiffness;
      elements.push_back(e);
      std::swap(e.rowIndex, e.columnIndex);
      e.value = linkData.derivative * particleInvMass[1] * -linkStiffness;
      elements.push_back(e);

      e.rowIndex = particleIndices[0];
      e.columnIndex = particleIndices[0];
      e.value = linkData.derivative * particleInvMass[0] * linkStiffness;
      elements.push_back(e);
      e.rowIndex = particleIndices[1];
      e.columnIndex = particleIndices[1];
      e.value = linkData.derivative * particleInvMass[1] * linkStiffness;
      elements.push_back(e);

      constraintDeltas[particleIndices[0]] -= linkData.delta * particleInvMass[0] * -linkStiffness;
      constraintDeltas[particleIndices[1]] += linkData.delta * particleInvMass[1] * -linkStiffness;
    }

    PhysicsData::ParticleTensorMatrix::SortElements(elements.data(), elements.size());
    constraintDerivatives.BuildFromSortedElements(elements.data(), elements.size(), PhysicsData::ParticleDimension{ particles.size() }, PhysicsData::ParticleDimension{ particles.size() });
  }
  void BuildParticleVelocities(ParticleGroup::Type particles, PhysicsData::DeltaMatrix& particleVelocities)
  {
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();

    particleVelocities.BuildEmpty(PhysicsData::ParticleDimension{ particles.size() }, PhysicsData::OneDimension());
    for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
    {
      particleVelocities.AppendRow(particleIndex);
#if !defined(POSITION_BASED)
      particleVelocities.AppendTerm(0, particleComponents[particleIndex].velocity);
#endif
    }
  }

  //try implicit position-based
  void SolveLinksImplicitMultigrid(ParticleGroup::Type particles, LinkGroup::Type links, float timeStep, legit::CpuProfiler& profiler)
  {
    PhysicsData::ParticleTensorMatrix constraintDerivatives;
    PhysicsData::ParticleTensorMatrix identity;
    PhysicsData::ParticleTensorMatrix systemMatrix;
    std::vector<glm::vec2> constraintDeltas;
    std::vector<glm::vec2> rightSide;
    std::vector<glm::vec2> derivativesTimesVelocity;
    PhysicsData::DeltaMatrix derivativesTimesVelocitySparse;
    PhysicsData::DeltaMatrix particleVelocities;
    constraintDeltas.resize(particles.size());
    derivativesTimesVelocity.resize(particles.size());
    rightSide.resize(particles.size());
    size_t iterationsCount = 10;

    {

      auto* particleComponents = particles.raw<ParticleComponent>();
      /* for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
       {
         particleComponents[particleIndex].pos += glm::vec2((nrand() - 0.5f) * 50.0f, (nrand() - 0.5f) * 50.0f);
       }*/

      {
        auto physicsTask = profiler.StartScopedTask("[Physics] Matrices", legit::Colors::greenSea);

        BuildConstraintMatrices(particles, links, constraintDerivatives, constraintDeltas.data());
      }

      /*{
        std::vector<glm::vec2> testDeltas;
        testDeltas.resize(particles.size());
        for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
        {
          testDeltas[particleIndex] = glm::vec2((nrand() - 0.5f) * 0.01f, (nrand() - 0.5f) * 0.01f);
        }

        std::vector<glm::vec2> predictedConstraintDeltas;
        predictedConstraintDeltas.resize(particles.size());
        BuildDenseVectorProduct<TensorVectorProduct>(constraintDerivatives, testDeltas.data(), predictedConstraintDeltas.data(), implicitStorage);

        for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
        {
          particleComponents[particleIndex].pos += testDeltas[particleIndex];
        }

        std::vector<glm::vec2> constraintDeltas2;
        constraintDeltas2.resize(particles.size());
        PhysicsData::ParticleTensorMatrix constraintDerivatives2;
        BuildConstraintMatrices(particles, links, constraintDerivatives2, constraintDeltas2.data());

        std::vector<glm::vec2> errs;
        errs.resize(particles.size());
        std::vector<glm::vec2> actualConstraintDeltas;
        actualConstraintDeltas.resize(particles.size());

        for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
        {
          actualConstraintDeltas[particleIndex] = constraintDeltas2[particleIndex] - constraintDeltas[particleIndex];
          errs[particleIndex] = actualConstraintDeltas[particleIndex] - predictedConstraintDeltas[particleIndex];
        }
        int p = 1;
      }*/

      {
        auto physicsTask = profiler.StartScopedTask("[Physics] System", legit::Colors::sunFlower);
        identity.BuildFromDiag(Tensor2f(1.0f / timeStep), PhysicsData::ParticleDimension{ particles.size() });
        systemMatrix.BuildFromSum(identity, 1.0f, constraintDerivatives, -timeStep, implicitStorage);

        BuildParticleVelocities(particles, particleVelocities);
        derivativesTimesVelocitySparse.BuildFromSparseProduct<TensorVectorProduct>(constraintDerivatives, particleVelocities, implicitStorage);

        BuildDenseVectorFromSparse(derivativesTimesVelocity.data(), derivativesTimesVelocitySparse);

        for (size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
        {
          rightSide[particleIndex] = constraintDeltas[particleIndex] + derivativesTimesVelocity[particleIndex] * timeStep;
        }

        if (1)
        {
          std::vector<glm::vec2> velocityDeltas;
          velocityDeltas.resize(particles.size());

          IterateGaussSeidel<TensorVectorProduct>(systemMatrix, rightSide.data(), velocityDeltas.data(), 10, implicitStorage);
          ApplyVelocityDenseDeltas(particles, velocityDeltas.data());
        }
        else
        {
          implicitMultigridSolver.LoadSystem(systemMatrix, implicitStorage);
          implicitMultigridSolver.Solve(rightSide.data(), 10, implicitStorage);
          ApplyVelocityDenseDeltas(particles, implicitMultigridSolver.GetValues());
        }
      }
    }

    {
      auto physicsTask = profiler.StartScopedTask("[Physics] Load", legit::Colors::pomegranate);
      //multigridSolver.LoadSystem(systemMatrix, links.size(), implicitStorage);
    }
  }

  struct Collision
  {
    size_t particleIndices[2];
    glm::vec2 norm;
    float depth;
  };

  void BuildSparseJacobian(ParticleGroup::Type particles, Collision* collisions, size_t collisionsCount,
    PhysicsData::JacobianMatrix& jacobianMatrix, float* positionRightSideVector, float* velocityRightSideVector, float* accelerationRightSide)
  {
    jacobianMatrix.BuildEmpty(PhysicsData::JointDimension{ collisionsCount }, PhysicsData::ParticleDimension{ particles.size() });
    ParticleComponent* particleComponents = particles.raw<ParticleComponent>();
    MassComponent* massComponents = particles.raw<MassComponent>();

    for (size_t collisionIndex = 0; collisionIndex < collisionsCount; collisionIndex++)
    {
      size_t particleIndex0 = collisions[collisionIndex].particleIndices[0];
      MassComponent& massComponent0 = massComponents[particleIndex0];
      ParticleComponent& particleComponent0 = particleComponents[particleIndex0];

      size_t particleIndex1 = collisions[collisionIndex].particleIndices[1];
      MassComponent& massComponent1 = massComponents[particleIndex1];
      ParticleComponent& particleComponent1 = particleComponents[particleIndex1];

      float deltaPos = collisions[collisionIndex].depth;

      float deltaAcceleration = -glm::dot(collisions[collisionIndex].norm, particleComponent0.acceleration - particleComponent1.acceleration);
#if defined(POSITION_BASED)
      float deltaVelocity = 0.0f;
#else
      float deltaVelocity = -glm::dot(collisions[collisionIndex].norm, particleComponent0.velocity - particleComponent1.velocity);
#endif

      {
        jacobianMatrix.AppendRow(collisionIndex);
        jacobianMatrix.AppendTerm(particleIndex0, collisions[collisionIndex].norm);
        jacobianMatrix.AppendTerm(particleIndex1, -collisions[collisionIndex].norm);
      }
      {
        positionRightSideVector[collisionIndex] = deltaPos;
      }
      {
        velocityRightSideVector[collisionIndex] = deltaVelocity;
      }
      {
        accelerationRightSide[collisionIndex] = deltaAcceleration;
      }
    }
    //jacobianMatrix.Finalize();
    assert(jacobianMatrix.CheckSortedIndices());
  }



  void SolveCollisions(ParticleGroup::Type particles, legit::CpuProfiler& profiler)
  {
    auto* particleComponents = particles.raw<ParticleComponent>();
    float radius0 = 5.0f;
    float radius1 = 5.0f;
    std::vector<Collision> collisions;
    collisions.clear();

    size_t groundParticleIndex = 0;

    for (size_t particleIndex0 = groundParticleIndex + 1; particleIndex0 < particles.size(); particleIndex0++)
    {
      glm::vec2 pos0 = particleComponents[particleIndex0].pos;
      float groundLevel = -300.0f;
      if (pos0.y - radius0 < groundLevel)
      {
        Collision collision;
        collision.particleIndices[0] = groundParticleIndex;
        collision.particleIndices[1] = particleIndex0;
        collision.norm = glm::vec2(0.0f, -1.0f);
        collision.depth = -(pos0.y - radius0) + groundLevel;
        collisions.push_back(collision);
      }
      for (size_t particleIndex1 = particleIndex0 + 1; particleIndex1 < particles.size(); particleIndex1++)
      {
        glm::vec2 pos1 = particleComponents[particleIndex1].pos;
        glm::vec2 delta = pos0 - pos1;
        float sqrDist = glm::dot(delta, delta);
        if (sqrDist < (radius0 + radius1) * (radius0 + radius1))
        {
          Collision collision;
          collision.particleIndices[0] = particleIndex0;
          collision.particleIndices[1] = particleIndex1;
          collision.norm = glm::normalize(delta);
          collision.depth = (radius0 + radius1) - sqrt(sqrDist);
          collisions.push_back(collision);
        }
      }
    }
    PhysicsData::JacobianMatrix jacobianMatrix;
    std::vector<float> positionRightSide;
    positionRightSide.resize(collisions.size());
    std::vector<float> velocityRightSide;
    velocityRightSide.resize(collisions.size());
    std::vector<float> accelerationRightSide;
    accelerationRightSide.resize(collisions.size());
    BuildSparseJacobian(particles, collisions.data(), collisions.size(), jacobianMatrix, positionRightSide.data(), velocityRightSide.data(), accelerationRightSide.data());

    PhysicsData::MassMatrix massMatrix;
    PhysicsData::JacobianTransposedMatrix jacobianTransposedMatrix;

    jacobianTransposedMatrix.BuildFromTransposed(jacobianMatrix, stackStorage);
    BuildSparseMassMatrix(particles, massMatrix);
    massJacobianTransposed.BuildFromSparseProduct<VectorScalarProduct>(massMatrix, jacobianTransposedMatrix, stackStorage);
    systemMatrix.BuildFromSparseProduct<DotProduct>(jacobianMatrix, massJacobianTransposed, stackStorage);
    {
      auto physicsTask = profiler.StartScopedTask("[Physics] Load", legit::Colors::pomegranate);
      multigridSolver.LoadSystem(systemMatrix, stackStorage);
    }
    {
      size_t iterationsCount = 5;
      auto positionLambdasSparseHandle = stackStorage.GetHandle<PhysicsData::RightSideMatrix>();
      auto positionDeltaMatrixHandle = stackStorage.GetHandle<PhysicsData::DeltaMatrix>();

      auto& positionLambdasSparse = positionLambdasSparseHandle.Get();
      auto& positionDeltaMatrix = positionDeltaMatrixHandle.Get();
      if (0)
      {
        auto positionLambdasHandle = stackStorage.GetHandle<std::vector<float>>();
        auto& positionLambdas = positionLambdasHandle.Get();
        positionLambdas.clear();
        positionLambdas.resize(collisions.size());

        IterateGaussSeidel<BasicSpace<float>>(systemMatrix, positionRightSide.data(), positionLambdas.data(), iterationsCount, implicitStorage);
        BuildSparseVectorFromDense(positionLambdas.data(), positionLambdas.size(), positionLambdasSparse);
      }
      else
      {
        multigridSolver.Solve(positionRightSide.data(), iterationsCount, stackStorage);
        BuildSparseVectorFromDense(multigridSolver.GetValues(), collisions.size(), positionLambdasSparse);
      }

      positionDeltaMatrix.BuildFromSparseProduct<VectorScalarProduct>(massJacobianTransposed, positionLambdasSparse, stackStorage);
      ApplyPositionDeltas(particles, positionDeltaMatrix);
    }
  }
}