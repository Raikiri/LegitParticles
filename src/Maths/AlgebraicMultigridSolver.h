#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "StackStorage.h"
#include "SparseMatrix.h"
#include "SparseMatrixUtils.h"
#include "TensorMaths.h"

namespace almost
{
  template<typename Space, typename DimensionType>
  struct MultigridLayer
  {
    using ValueType0 = typename Space::ValueType0;
    using ValueType1 = typename Space::ValueType1;
    using ResultType = typename Space::ResultType;
    using ScalarType = typename Space::ScalarType;
    using FineDimensionType = typename DimensionType;
    SparseMatrix<DimensionType, DimensionType, ValueType0> systemMatrix;

    size_t valuesCount;
    std::vector<ResultType> rightSide;
    std::vector<ValueType1> values;
    std::vector<ResultType> residuals;

    SparseMatrix<DimensionType, FineDimensionType, ValueType0> tmpMatrix;

    SparseMatrix<FineDimensionType, DimensionType, ScalarType> restrictionMatrix; //fine prev->this
    SparseMatrix<DimensionType, FineDimensionType, ScalarType> interpolationMatrix; //this->fine prev

    void Resize(size_t _valuesCount)
    {
      this->valuesCount = _valuesCount;
      rightSide.resize(valuesCount);
      values.resize(valuesCount);
      residuals.resize(valuesCount);
    }
    template<typename StorageType>
    void IterateGaussSeidel(size_t iterationsCount, StorageType &storage)
    {
      almost::IterateGaussSeidel<Space>(systemMatrix, rightSide.data(), values.data(), iterationsCount, storage);
    }

    template<typename StorageType>
    void BuildResiduals(StorageType &storage)
    {
      almost::BuildResiduals<Space>(systemMatrix, rightSide.data(), values.data(), residuals.data(), storage);
    }

    void ResetValues()
    {
      std::fill_n(values.data(), values.size(), ValueType1(0));
    }


    /*template<typename StorageType>
    void Downsample(const MultigridLayer &fineSolver, StorageType & storage)
    {
      this->systemMatrix = this->restrictionMatrix.SparseMatrixProduct<ValueType>(fineSolver.systemMatrix.SparseMatrixProduct<ValueType>(fineSolver.interpolationMatrix));

      SparseMatrix<FineIndexType, char, ValueType> fineResidualsSparse;
      fineSolver.DenseVectorToSparseMatrix(fineSolver.residuals.data(), fineResidualsSparse);
      
      rightSide.resize(valuesCount);
      SparseMatrix<IndexType, char, ValueType> rightSideSparse = this->restrictionMatrix.SparseMatrixProduct<ValueType>(fineResidualsSparse);
      //assert(rightSideSparse.GetRowsCount() == rightSide.size());
      this->SparseMatrixToDenseVector(rightSideSparse, rightSide.data());
    }

    template<typename StorageType>
    void Upsample(const MultigridLayer &coarseSolver, StorageType& storage)
    {
      SparseMatrix<CoarseIndexType, char, ValueType> coarseValuesSparse;
      coarseSolver.DenseVectorToSparseMatrix(coarseSolver.values.data(), coarseValuesSparse);

      SparseMatrix<IndexType, char, ValueType> valuesSparse;
      
      valuesSparse.BuildFromSparseMatrixProduct(interpolationMatrix, coarseValuesSparse, storage);
      std::vector<ValueType> valuesDelta;
      valuesDelta.resize(values.size());
      SparseMatrixToDenseVector(valuesSparse, valuesDelta.data());
      for (size_t valueIndex = 0; valueIndex < values.size(); valueIndex++)
        values[valueIndex] += valuesDelta[valueIndex];
    }*/
  };

  template<typename Space, typename DimensionType>
  struct AlgebraicMultigridSolver
  {
    using ValueType0 = typename Space::ValueType0;
    using ValueType1 = typename Space::ValueType1;
    using ResultType = typename Space::ResultType;

    AlgebraicMultigridSolver()
    {
      layerSolvers.resize(4);
    }

    template<typename StorageType>
    void LoadSystem(const SparseMatrix<DimensionType, DimensionType, ValueType0>& systemMatrix, StorageType &storage)
    {
      assert(systemMatrix.rowDimension.size == systemMatrix.columnDimension.size);
      layerSolvers[0].Resize(systemMatrix.rowDimension.size);
      layerSolvers[0].systemMatrix = systemMatrix;

      for (size_t i = 0; i < layerSolvers.size() - 1; i++)
      {
        almost::BuildInterpolationMatrix<Space>(layerSolvers[i].systemMatrix, layerSolvers[i].valuesCount, layerSolvers[i].interpolationMatrix, storage);
        layerSolvers[i].restrictionMatrix.BuildFromTransposed(layerSolvers[i].interpolationMatrix, storage);

        layerSolvers[i + 1].Resize(layerSolvers[i].restrictionMatrix.rowDimension.size);
        layerSolvers[i].tmpMatrix.BuildFromSparseProduct<Space>(layerSolvers[i].restrictionMatrix, layerSolvers[i].systemMatrix, storage);
        layerSolvers[i + 1].systemMatrix.BuildFromSparseProduct<Space>(layerSolvers[i].tmpMatrix, layerSolvers[i].interpolationMatrix, storage);
      }
    }

    template<typename StorageType>
    void Solve(const ResultType *rightSide, size_t iterationsCount, StorageType &storage)
    {
      std::copy_n(rightSide, layerSolvers[0].valuesCount, layerSolvers[0].rightSide.data());
      for (size_t i = 1; i < layerSolvers.size(); i++)
      {
        layerSolvers[i - 1].ResetValues();
        layerSolvers[i - 1].IterateGaussSeidel(iterationsCount / 2, storage);
        layerSolvers[i - 1].BuildResiduals(storage);
        BuildDenseVectorProduct<Space>(layerSolvers[i - 1].restrictionMatrix, layerSolvers[i - 1].residuals.data(), layerSolvers[i].rightSide.data(), storage);
      }
      layerSolvers.back().ResetValues();

      int p = 1;
      for (size_t i = 0; i < layerSolvers.size(); i++)
      {
        size_t layerIndex = layerSolvers.size() - 1 - i;
        layerSolvers[layerIndex].IterateGaussSeidel(iterationsCount / 2, storage);

        if (layerIndex > 0)
        {
          std::vector<ValueType1> valuesDelta;
          valuesDelta.resize(layerSolvers[layerIndex - 1].valuesCount);

          almost::BuildDenseVectorProduct<Space>(layerSolvers[layerIndex - 1].interpolationMatrix, layerSolvers[layerIndex].values.data(), valuesDelta.data(), storage);
          for (size_t valueIndex = 0; valueIndex < layerSolvers[layerIndex - 1].valuesCount; valueIndex++)
            layerSolvers[layerIndex - 1].values[valueIndex] += valuesDelta[valueIndex];
        }
      }
    }

    ValueType1 *GetValues()
    {
      return layerSolvers[0].values.data();
    }

  private:
    std::vector<MultigridLayer<Space, DimensionType>> layerSolvers;
  };
}