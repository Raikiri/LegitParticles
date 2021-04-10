#pragma once
#include "SparseMatrix.h"
namespace almost
{
  template<typename DimensionType, typename ValueType>
  static void ComputeSparseDiag(const SparseMatrix<DimensionType, DimensionType, ValueType>& matrix, ValueType* diagCoeffs)
  {
    using IndexType = typename DimensionType::IndexType;
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

  template<typename Space, typename DimensionType, typename StorageType>
  static void BuildSparseConnectivityMatrix(const SparseMatrix<DimensionType, DimensionType, typename Space::ValueType0>& srcMatrix, typename Space::ValueType0* srcDiag, SparseMatrix<DimensionType, DimensionType, typename Space::ScalarType>& dstConnectivityMatrix, StorageType& tmpStorage)
  {
    using IndexType = typename DimensionType::IndexType;
    //BuildSparseDiag(srcMatrix, dstDiagCoeff);
    dstConnectivityMatrix.BuildEmpty(srcMatrix.rowDimension, srcMatrix.columnDimension);

    auto adjacentRadiiHandle = tmpStorage.template GetHandle < std::vector<typename Space::ScalarType> >();
    auto adjacentRadii = adjacentRadiiHandle.Get();

    for (size_t valueNumber = 0; valueNumber < srcMatrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = srcMatrix.GetRowIndex(valueNumber);

      size_t adjacentValuesCount = srcMatrix.GetRowTermsCount(valueNumber);
      adjacentRadii.resize(adjacentValuesCount);
      const auto* rowTerms = srcMatrix.GetRowTerms(valueNumber);
      float maxRadius = 0.0f;
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
        //assert(srcDiag[rowIndex] > 0.0f);
        //assert(srcDiag[columnIndex] > 0.0f);
        float radius = Space::SpectralRadius(rowTerms[adjacentValueNumber].value, srcDiag[rowIndex], srcDiag[columnIndex]);
        adjacentRadii[adjacentValueNumber] = radius;
        if (columnIndex != rowIndex)
          maxRadius = std::max(maxRadius, radius);
      }

      dstConnectivityMatrix.AppendRow(rowIndex);
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
        if (adjacentRadii[adjacentValueNumber] > maxRadius * 0.5 || columnIndex == rowIndex)
        {
          dstConnectivityMatrix.AppendTerm(rowTerms[adjacentValueNumber].columnIndex, 1.0f);
        }
      }
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


  template<typename FineDimensionType, typename CoarseDimensionType, typename ScalarType, typename StorageType>
  static void BuildTentativeInterpolationMatrix(
    const SparseMatrix<FineDimensionType, FineDimensionType, ScalarType>& srcFineConnectivityMatrix,
    typename FineDimensionType::IndexType valuesCount,
    SparseMatrix<CoarseDimensionType, FineDimensionType, ScalarType>& dstInterpolationMatrix,
    StorageType& storage)
  {
    using CoarseIndexType = typename CoarseDimensionType::IndexType;
    using FineIndexType = typename FineDimensionType::IndexType;
    auto aggregationIndicesHandle = storage.template GetHandle<std::vector<CoarseIndexType>>();
    auto aggregationIndices = aggregationIndicesHandle.Get();

    aggregationIndices.resize(valuesCount);

    std::fill(aggregationIndices.begin(), aggregationIndices.end(), CoarseIndexType(-1));
    CoarseIndexType aggregatesCount = 0;
    for (size_t valueNumber = 0; valueNumber < srcFineConnectivityMatrix.GetRowsCount(); valueNumber++)
    {
      FineIndexType valueIndex = srcFineConnectivityMatrix.GetRowIndex(valueNumber);
      size_t adjacentValuesCount = srcFineConnectivityMatrix.GetRowTermsCount(valueNumber);
      const auto *rowTerms = srcFineConnectivityMatrix.GetRowTerms(valueNumber);
      bool isUnassigned = true;
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        FineIndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
        if (aggregationIndices[columnIndex] != CoarseIndexType(-1))
          isUnassigned = false;
      }
      if (isUnassigned)
      {
        for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
        {
          FineIndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
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
      const auto *rowTerms = srcFineConnectivityMatrix.GetRowTerms(valueNumber);
      CoarseIndexType dstAggregationIndex = CoarseIndexType(-1);
      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        FineIndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
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

    dstInterpolationMatrix.BuildEmpty(srcFineConnectivityMatrix.rowDimension, CoarseDimensionType{ aggregatesCount });
    //SparseMatrix<IndexType, IndexType, ValueType> tentativeInterpolationMatrix;
    for (size_t valueIndex = 0; valueIndex < valuesCount; valueIndex++)
    {
      FineIndexType rowIndex = valueIndex;
      CoarseIndexType columnIndex = aggregationIndices[valueIndex];
      assert(columnIndex != CoarseIndexType(-1));
      ScalarType coeff = ScalarType(1) / sqrt(ScalarType(aggregateSizes[columnIndex]));
      dstInterpolationMatrix.AppendRow(rowIndex);
      dstInterpolationMatrix.AppendTerm(columnIndex, coeff);
    }
    assert(dstInterpolationMatrix.CheckSortedIndices());
  }

  template<typename Space, typename FineDimensionType, typename CoarseDimensionType, typename StorageType>
  static void BuildInterpolationMatrix(const SparseMatrix<FineDimensionType, FineDimensionType, typename Space::ValueType0> &srcFineSystemMatrix, typename FineDimensionType::IndexType valuesCount, SparseMatrix<FineDimensionType, CoarseDimensionType, typename Space::ScalarType> &dstInterpolationMatrix, StorageType &storage)
  {
    using ConnectivityMatrix = SparseMatrix<FineDimensionType, FineDimensionType, typename Space::ScalarType>;
    auto diagCoeffsHandle = storage.template GetHandle<std::vector<typename Space::ValueType0>>();
    auto diagCoeffs = diagCoeffsHandle.Get();
    diagCoeffs.resize(valuesCount);

    auto connectivityMatrixHandle = storage.template GetHandle< ConnectivityMatrix >();
    auto connectivityMatrix = connectivityMatrixHandle.Get();

    almost::ComputeSparseDiag(srcFineSystemMatrix, diagCoeffs.data());
    almost::BuildSparseConnectivityMatrix<Space>(srcFineSystemMatrix, diagCoeffs.data(), connectivityMatrix, storage);


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

  template<typename Space, typename MatrixType, typename ValueType, typename StorageType>
  void IterateGaussSeidel(const MatrixType systemMatrix, const ValueType *rightSide, ValueType *values, size_t iterationsCount, StorageType &storage)
  {
    using IndexType = typename MatrixType::RowIndexType;
    for (size_t iterationIndex = 0; iterationIndex < iterationsCount; iterationIndex++)
    {
      for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
      {
        IndexType rowIndex = systemMatrix.GetRowIndex(valueNumber);
        size_t valueIndex = rowIndex;
        assert(valueNumber == rowIndex); //typically we have non-zero at least on diag

        size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);

        const auto *rowTerms = systemMatrix.GetRowTerms(valueNumber);
        ValueType remainder = rightSide[valueIndex];
        typename MatrixType::ValueType ownCoeff = typename MatrixType::ValueType(-1.0f);

        for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
        {
          IndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
          size_t adjacentValueIndex = columnIndex;

          if (columnIndex == rowIndex)
            ownCoeff = rowTerms[adjacentValueNumber].value;
          else
            remainder -= Space::Mul(rowTerms[adjacentValueNumber].value, values[adjacentValueIndex]);
        }
        //assert(ownCoeff > 0.0f);
        values[valueIndex] = Space::Divide(remainder, ownCoeff);
      }
    }
  }

  template<typename Space, typename DimensionType, typename StorageType>
  void BuildResiduals(
    const SparseMatrix<DimensionType, DimensionType, typename Space::ValueType0> systemMatrix,
    const typename Space::ResultType* rightSide,
    const typename Space::ValueType1* values,
    typename Space::ResultType *residuals,
    StorageType& storage)
  {
    using IndexType = typename DimensionType::IndexType;

    for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
    {
      IndexType rowIndex = systemMatrix.GetRowIndex(valueNumber);
      size_t valueIndex = rowIndex; //valueNumber

      size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);
      const auto *rowTerms = systemMatrix.GetRowTerms(valueNumber);
      typename Space::ResultType remainder = rightSide[valueIndex];

      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        IndexType columnIndex = rowTerms[adjacentValueNumber].columnIndex;
        size_t adjacentValueIndex = columnIndex;
        remainder -= Space::Mul(rowTerms[adjacentValueNumber].value, values[adjacentValueIndex]);
      }
      residuals[valueIndex] = remainder;
    }
  }

  template<typename Space, typename MatrixType, typename ValueType, typename StorageType>
  void BuildDenseVectorProduct(const MatrixType systemMatrix, const ValueType* denseVector, ValueType* denseProduct, StorageType& storage)
  {
    using ColumnIndexType = typename MatrixType::ColumnIndexType;

    for (size_t valueNumber = 0; valueNumber < systemMatrix.GetRowsCount(); valueNumber++)
    {
      size_t rowIndex = systemMatrix.GetRowIndex(valueNumber);
      size_t valueIndex = rowIndex; //valueNumber

      size_t adjacentValuesCount = systemMatrix.GetRowTermsCount(valueNumber);
      const auto *rowTerms = systemMatrix.GetRowTerms(valueNumber);
      ValueType sum = ValueType(0);

      for (size_t adjacentValueNumber = 0; adjacentValueNumber < adjacentValuesCount; adjacentValueNumber++)
      {
        size_t columnIndex = rowTerms[adjacentValueNumber].columnIndex;
        size_t adjacentValueIndex = columnIndex;
        sum += Space::Mul(rowTerms[adjacentValueNumber].value, denseVector[adjacentValueIndex]);
      }
      denseProduct[valueIndex] = sum;
    }
  }
}