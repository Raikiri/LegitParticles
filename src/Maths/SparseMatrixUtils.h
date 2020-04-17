#pragma once
#include "SparseMatrix.h"
namespace almost
{
  template<typename IndexType, typename ValueType>
  static void ComputeSparseDiag(const SparseMatrix<IndexType, IndexType, ValueType>& matrix, ValueType* diagCoeffs)
  {
    //diagCoeffs.resize(valuesCount);

    for (size_t valueNumber = 0; valueNumber < matrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = matrix.GetRowIndex(valueNumber);
      IndexType columnIndex = rowIndex;
      ValueType diagCoeff;
      if (matrix.FindTerm(valueNumber, columnIndex, diagCoeff))
        diagCoeffs[rowIndex] = diagCoeff;
      else
        diagCoeffs[rowIndex] = ValueType(0);
    }
  }

  template<typename IndexType, typename ValueType, typename StorageType>
  static void BuildSparseConnectivityMatrix(const SparseMatrix<IndexType, IndexType, ValueType>& srcMatrix, ValueType* srcDiag, SparseMatrix<IndexType, IndexType, ValueType>& dstConnectivityMatrix, StorageType& tmpStorage)
  {
    //BuildSparseDiag(srcMatrix, dstDiagCoeff);
    dstConnectivityMatrix.Clear();

    //std::vector<IndexType> connectedIndices;
    //std::vector<float> connectedWeights;
    auto adjacentIndicesHandle = tmpStorage.template GetHandle < std::vector<IndexType> >();
    auto adjacentIndices = adjacentIndicesHandle.Get();
    auto adjacentWeightsHandle = tmpStorage.template GetHandle < std::vector<ValueType> >();
    auto adjacentWeights = adjacentWeightsHandle.Get();
    auto adjacentRadiiHandle = tmpStorage.template GetHandle < std::vector<ValueType> >();
    auto adjacentRadii = adjacentRadiiHandle.Get();

    auto connectedIndicesHandle = tmpStorage.template GetHandle < std::vector<IndexType> >();
    auto connectedIndices = connectedIndicesHandle.Get();
    auto connectedWeightsHandle = tmpStorage.template GetHandle < std::vector<ValueType> >();
    auto connectedWeights = connectedWeightsHandle.Get();

    for (size_t valueNumber = 0; valueNumber < srcMatrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = srcMatrix.GetRowIndex(valueNumber);

      size_t adjacentValuesCount = srcMatrix.GetRowTermsCount(valueNumber);
      adjacentIndices.resize(adjacentValuesCount);
      adjacentWeights.resize(adjacentValuesCount);
      adjacentRadii.resize(adjacentValuesCount);
      srcMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());
      ValueType maxRadius = 0.0f;
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = adjacentIndices[adjacentValueNumber];
        ValueType coeff = abs(adjacentWeights[adjacentValueNumber]);
        assert(srcDiag[rowIndex] > 0.0f);
        assert(srcDiag[columnIndex] > 0.0f);
        ValueType radius = coeff / (sqrt(srcDiag[rowIndex]) * sqrt(srcDiag[columnIndex]));
        adjacentRadii[adjacentValueNumber] = radius;
        if (columnIndex != rowIndex)
          maxRadius = std::max(maxRadius, radius);
      }

      connectedWeights.clear();
      connectedIndices.clear();
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        if (adjacentRadii[adjacentValueNumber] > maxRadius * 0.5f)
        {
          connectedWeights.push_back(1.0f);
          connectedIndices.push_back(adjacentIndices[adjacentValueNumber]);
        }
      }
      dstConnectivityMatrix.AddRow(rowIndex, connectedIndices.data(), connectedWeights.data(), connectedIndices.size());
    }
    assert(dstConnectivityMatrix.CheckSortedIndices());
    //interpolationMatrix = coarseSolver.restrictionMatrix.GetTransposed();

    /*struct Values {};
    using ValuesDim = almost::StaticDimension<Values, 3>;
    struct NullspaceVectors {};
    using NullspaceVectorsDim = almost::StaticDimension<NullspaceVectors, 7>;
    struct Basis {};
    using BasisDim = almost::DynamicDimension<Basis>;

    ValuesDim valuesDim;
    NullspaceVectorsDim nullspaceDim;
    BasisDim basisDim(1); //size specified in runtime


    auto valueToNullspace = almost::MakeStaticTensor<float, ValuesDim, NullspaceVectorsDim>(); //on stack, static array
    auto nullspaceToBasis = almost::MakeDynamicTensor<float>(nullspaceDim, basisDim); //basisDim is dynamic, allocate on heap
    auto valuesInBasis = almost::MakeDynamicTensor<float>(valuesDim, basisDim); //basisDim is dynamic, allocate on heap
    for (auto i : valuesDim)
      for (auto k : basisDim)
        valuesInBasis.Get(i, k) = (i.value == k.value) ? 1.0f : -42.0f;

    for (auto i : valuesDim)
      for (auto j : nullspaceDim)
        for (auto k : basisDim)
          valuesInBasis.Get(i, k) += valueToNullspace.Get(i, j) * nullspaceToBasis.Get(j, k);*/


          /*for (auto i : valuesDim)
          {
            for (auto k : basisDim)
            {
              valuesInBasis.Get(i, k) = (i.value == k.value) ? 1.0f : -42.0f;
            }
          }*/

          /*StaticMatrix<3, 3, float> test;
          test.Get(0, 0) = 12; test.Get(0, 1) = -51; test.Get(0, 2) = 4;
          test.Get(1, 0) = 6; test.Get(1, 1) =  167; test.Get(1, 2) = -68;
          test.Get(2, 0) =-4; test.Get(2, 1) = 24; test.Get(2, 2) = -41;

          StaticMatrix<3, 3, float> q;
          StaticMatrix<3, 3, float> r;
          MatrixDecomposeQR(test, q, r);*/
    int pp = 1;
  }
  template<typename IndexType, typename ValueType>
  static void DenseVectorToSparseMatrix(const ValueType* denseVector, size_t valuesCount, SparseMatrix<IndexType, char, ValueType>& dstSparseMatrix)
  {
    dstSparseMatrix.Clear();
    size_t rowsCount = valuesCount;
    for (size_t rowIndex = 0; rowIndex < rowsCount; rowIndex++)
    {
      char columnIndex = 0;
      float value = denseVector[rowIndex];
      if (value != 0)
        dstSparseMatrix.AddRow(rowIndex, &columnIndex, &value, 1);
    }
  }
  template<typename IndexType, typename ValueType>
  static void SparseMatrixToDenseVector(const SparseMatrix<IndexType, char, ValueType>& sparseMatrix, ValueType* denseVector)
  {
    size_t rowsCount = sparseMatrix.GetRowsCount();
    //assert(sparseMatrix.GetRowsCount() == valuesCount);
    for (size_t rowNumber = 0; rowNumber < rowsCount; rowNumber++)
    {
      assert(sparseMatrix.GetRowTermsCount(rowNumber) == 1);
      size_t rowIndex = sparseMatrix.GetRowIndex(rowNumber);

      ValueType value;
      char columnIndex;
      sparseMatrix.GetRowTerms(rowNumber, &columnIndex, &value);
      assert(columnIndex == 0);
      denseVector[rowIndex] = value;
    }
  }


  template<typename FineIndexType, typename CoarseIndexType, typename ValueType, typename StorageType>
  static void BuildTentativeInterpolationMatrix(const SparseMatrix<FineIndexType, FineIndexType, ValueType>& srcFineConnectivityMatrix, FineIndexType valuesCount, SparseMatrix<CoarseIndexType, FineIndexType, ValueType>& dstInterpolationMatrix, StorageType& storage)
  {
    auto aggregationIndicesHandle = storage.template GetHandle<std::vector<CoarseIndexType>>();
    auto aggregationIndices = aggregationIndicesHandle.Get();
    auto adjacentIndicesHandle = storage.template GetHandle<std::vector<FineIndexType>>();
    auto adjacentIndices = adjacentIndicesHandle.Get();
    auto adjacentWeightsHandle = storage.template GetHandle<std::vector<ValueType>>();
    auto adjacentWeights = adjacentWeightsHandle.Get();

    aggregationIndices.resize(valuesCount);

    std::fill(aggregationIndices.begin(), aggregationIndices.end(), CoarseIndexType(-1));
    CoarseIndexType aggregatesCount = 0;
    for (size_t valueNumber = 0; valueNumber < srcFineConnectivityMatrix.GetRowsCount(); valueNumber++)
    {
      FineIndexType valueIndex = srcFineConnectivityMatrix.GetRowIndex(valueNumber);
      size_t adjacentValuesCount = srcFineConnectivityMatrix.GetRowTermsCount(valueNumber);
      adjacentIndices.resize(adjacentValuesCount);
      adjacentWeights.resize(adjacentValuesCount);
      srcFineConnectivityMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());
      bool isUnassigned = true;
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        FineIndexType columnIndex = adjacentIndices[adjacentValueNumber];
        if (aggregationIndices[columnIndex] != CoarseIndexType(-1))
          isUnassigned = false;
      }
      if (isUnassigned)
      {
        for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
        {
          FineIndexType columnIndex = adjacentIndices[adjacentValueNumber];
          aggregationIndices[columnIndex] = aggregatesCount;
        }
        aggregatesCount++;
      }
    }

    for (size_t valueNumber = 0; valueNumber < srcFineConnectivityMatrix.GetRowsCount(); valueNumber++)
    {
      FineIndexType valueIndex = srcFineConnectivityMatrix.GetRowIndex(valueNumber);
      if (aggregationIndices[valueIndex] != CoarseIndexType(-1))
        continue;

      size_t adjacentValuesCount = srcFineConnectivityMatrix.GetRowTermsCount(valueNumber);
      adjacentIndices.resize(adjacentValuesCount);
      adjacentWeights.resize(adjacentValuesCount);
      srcFineConnectivityMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());
      CoarseIndexType dstAggregationIndex = CoarseIndexType(-1);
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        FineIndexType columnIndex = adjacentIndices[adjacentValueNumber];
        if (aggregationIndices[columnIndex] != CoarseIndexType(-1))
          dstAggregationIndex = aggregationIndices[columnIndex];
      }
      assert(dstAggregationIndex != CoarseIndexType(-1));
      aggregationIndices[valueIndex] = dstAggregationIndex;
    }

    auto aggregateSizesHandle = storage.template GetHandle< std::vector<size_t> >();
    auto aggregateSizes = aggregateSizesHandle.Get();
    aggregateSizes.resize(aggregatesCount);
    std::fill(aggregateSizes.begin(), aggregateSizes.end(), 0);
    for (auto aggregateIndex : aggregationIndices)
    {
      if (aggregateIndex != CoarseIndexType(-1))
        aggregateSizes[aggregateIndex]++;
    }

    dstInterpolationMatrix.Clear();
    //SparseMatrix<IndexType, IndexType, ValueType> tentativeInterpolationMatrix;
    for (size_t valueIndex = 0; valueIndex < valuesCount; valueIndex++)
    {
      FineIndexType rowIndex = valueIndex;
      CoarseIndexType columnIndex = aggregationIndices[valueIndex];
      assert(columnIndex != CoarseIndexType(-1));
      ValueType coeff = ValueType(1) / sqrt(ValueType(aggregateSizes[columnIndex]));
      dstInterpolationMatrix.AddRow(rowIndex, &columnIndex, &coeff, 1);
    }
    assert(dstInterpolationMatrix.CheckSortedIndices());
  }

  template<typename FineIndexType, typename CoarseIndexType, typename ValueType, typename StorageType>
  static void BuildInterpolationMatrix(const SparseMatrix<FineIndexType, FineIndexType, ValueType> &srcFineSystemMatrix, FineIndexType valuesCount, SparseMatrix<CoarseIndexType, FineIndexType, ValueType> &dstInterpolationMatrix, StorageType &storage)
  {
    using SystemMatrix = SparseMatrix<FineIndexType, FineIndexType, ValueType>;
    auto diagCoeffsHandle = storage.template GetHandle<std::vector<ValueType>>();
    auto diagCoeffs = diagCoeffsHandle.Get();
    diagCoeffs.resize(valuesCount);

    auto connectivityMatrixHandle = storage.template GetHandle< SystemMatrix>();
    auto connectivityMatrix = connectivityMatrixHandle.Get();

    almost::ComputeSparseDiag(srcFineSystemMatrix, diagCoeffs.data());
    almost::BuildSparseConnectivityMatrix(srcFineSystemMatrix, diagCoeffs.data(), connectivityMatrix, storage);

    dstInterpolationMatrix.Clear();

    BuildTentativeInterpolationMatrix(connectivityMatrix, valuesCount, dstInterpolationMatrix, storage);
    /*SparseMatrix<IndexType, IndexType, ValueType> smootherMatrix;
    for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = systemMatrix.GetRowIndex(valueNumber);

      size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);
      adjacentIndices.resize(adjacentValuesCount);
      adjacentWeights.resize(adjacentValuesCount);
      systemMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());

      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = adjacentIndices[adjacentValueNumber];
        float mult = 0.1f;
        float base = 0.0f;
        if (rowIndex == columnIndex)
        {
          base = 1.0f;
        }
        adjacentWeights[adjacentValueNumber] = base - mult * adjacentWeights[adjacentValueNumber] / diagCoeffs[rowIndex];
        //adjacentWeights[adjacentValueNumber] = base - mult * adjacentWeights[adjacentValueNumber];
      }
      smootherMatrix.AddRow(rowIndex, adjacentIndices.data(), adjacentWeights.data(), adjacentValuesCount);
    }
    //interpolationMatrix = smootherMatrix.SparseMatrixProduct<ValueType>(tentativeInterpolationMatrix);
    interpolationMatrix = tentativeInterpolationMatrix;*/

    /*auto testMatrix = interpolationMatrix;
    size_t newDimsCount = systemMatrix.GetRowsCount() / 2;
    for (int rowIndex = 0; rowIndex < int(newDimsCount); rowIndex++)
    {
      std::vector<IndexType> columnIndices;
      std::vector<float> columnWeights;
      for (size_t adjacentNumber = 0; adjacentNumber < 3; adjacentNumber++)
      {
        float weights[] = { 1, 2, 1 };
        int offsetIndex = int(rowIndex * 2) - 1 + int(adjacentNumber);
        if (offsetIndex >= 0 && offsetIndex < systemMatrix.GetRowsCount())
        {
          columnWeights.push_back(weights[adjacentNumber]);
          columnIndices.push_back(offsetIndex);
        }
      }
      testMatrix.AddRow(rowIndex, columnIndices.data(), columnWeights.data(), columnIndices.size());
    }*/
    //interpolationMatrix = testMatrix.GetTransposed();

    //auto transposedTest = tentativeInterpolationMatrix.GetTransposed();
    //auto test = transposedTest.SparseMatrixProduct<ValueType>(tentativeInterpolationMatrix);
    //int p = 1;
  }

  template<typename IndexType, typename ValueType, typename StorageType>
  void IterateGaussSeidel(const SparseMatrix<IndexType, IndexType, ValueType> systemMatrix, const ValueType *rightSide, ValueType *values, size_t iterationsCount, StorageType &storage)
  {
    auto adjacentIndicesHandle = storage.template GetHandle<std::vector<IndexType>>();
    auto adjacentIndices = adjacentIndicesHandle.Get();
    auto adjacentWeightsHandle = storage.template GetHandle<std::vector<ValueType>>();
    auto adjacentWeights = adjacentWeightsHandle.Get();

    for (size_t iterationIndex = 0; iterationIndex < iterationsCount; iterationIndex++)
    {
      for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
      {
        IndexType rowIndex = systemMatrix.GetRowIndex(valueNumber);
        size_t valueIndex = rowIndex;
        assert(valueNumber == rowIndex); //typically we have non-zero at least on diag

        size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);
        adjacentIndices.resize(adjacentValuesCount);
        adjacentWeights.resize(adjacentValuesCount);

        systemMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());
        float remainder = rightSide[valueIndex];
        float ownCoeff = -1.0f;

        for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
        {
          IndexType columnIndex = adjacentIndices[adjacentValueNumber];
          size_t adjacentValueIndex = columnIndex;

          if (columnIndex == rowIndex)
            ownCoeff = adjacentWeights[adjacentValueNumber];
          else
            remainder -= adjacentWeights[adjacentValueNumber] * values[adjacentValueIndex];
        }
        assert(ownCoeff > 0.0f);
        values[valueIndex] = remainder / ownCoeff;
      }
    }
  }

  template<typename IndexType, typename ValueType, typename StorageType>
  void BuildResiduals(const SparseMatrix<IndexType, IndexType, ValueType> systemMatrix, const ValueType* rightSide, const ValueType* values, ValueType *residuals, StorageType& storage)
  {
    auto adjacentIndicesHandle = storage.template GetHandle<std::vector<IndexType>>();
    auto adjacentIndices = adjacentIndicesHandle.Get();
    auto adjacentWeightsHandle = storage.template GetHandle<std::vector<ValueType>>();
    auto adjacentWeights = adjacentWeightsHandle.Get();

    for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = systemMatrix.GetRowIndex(valueNumber);
      size_t valueIndex = rowIndex; //valueNumber

      size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);
      adjacentIndices.resize(adjacentValuesCount);
      adjacentWeights.resize(adjacentValuesCount);

      systemMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());
      float remainder = rightSide[valueIndex];

      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = adjacentIndices[adjacentValueNumber];
        size_t adjacentValueIndex = columnIndex;
        remainder -= adjacentWeights[adjacentValueNumber] * values[adjacentValueIndex];
      }
      residuals[valueIndex] = remainder;
    }
  }

  template<typename IndexType, typename ValueType, typename StorageType>
  void BuildDenseVectorProduct(const SparseMatrix<IndexType, IndexType, ValueType> systemMatrix, const ValueType* denseVector, ValueType* denseProduct, StorageType& storage)
  {
    auto adjacentIndicesHandle = storage.template GetHandle<std::vector<IndexType>>();
    auto adjacentIndices = adjacentIndicesHandle.Get();
    auto adjacentWeightsHandle = storage.template GetHandle<std::vector<ValueType>>();
    auto adjacentWeights = adjacentWeightsHandle.Get();

    for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = systemMatrix.GetRowIndex(valueNumber);
      size_t valueIndex = rowIndex; //valueNumber

      size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);
      adjacentIndices.resize(adjacentValuesCount);
      adjacentWeights.resize(adjacentValuesCount);

      systemMatrix.GetRowTerms(valueNumber, adjacentIndices.data(), adjacentWeights.data());
      ValueType sum = ValueType(0);

      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = adjacentIndices[adjacentValueNumber];
        size_t adjacentValueIndex = columnIndex;
        sum += adjacentWeights[adjacentValueNumber] * denseVector[adjacentValueIndex];
      }
      denseProduct[valueIndex] = sum;
    }
  }
}