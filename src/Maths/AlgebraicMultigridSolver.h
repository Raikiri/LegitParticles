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
  template<typename IndexType, typename ValueType>
  struct MultigridLayer
  {
    using FineIndexType = IndexType;
    SparseMatrix<IndexType, IndexType, ValueType> systemMatrix;

    size_t valuesCount;
    std::vector<ValueType> rightSide;
    std::vector<ValueType> values;
    std::vector<ValueType> residuals;

    SparseMatrix<IndexType, FineIndexType, ValueType> tmpMatrix;

    SparseMatrix<FineIndexType, IndexType, ValueType> restrictionMatrix; //fine prev->this
    SparseMatrix<IndexType, FineIndexType, ValueType> interpolationMatrix; //this->fine prev

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
      almost::IterateGaussSeidel(systemMatrix, rightSide.data(), values.data(), iterationsCount, storage);
    }

    template<typename StorageType>
    void BuildResiduals(StorageType &storage)
    {
      almost::BuildResiduals(systemMatrix, rightSide.data(), values.data(), residuals.data(), storage);
    }

    void ResetValues()
    {
      std::fill_n(values.data(), values.size(), ValueType(0));
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

  template<typename IndexType, typename ValueType>
  struct AlgebraicMultigridSolver
  {
    AlgebraicMultigridSolver()
    {
      layerSolvers.resize(4);
    }

    template<typename StorageType>
    void LoadSystem(const SparseMatrix<IndexType, IndexType, ValueType>& systemMatrix, size_t valuesCount, StorageType &storage)
    {
      layerSolvers[0].Resize(valuesCount);
      layerSolvers[0].systemMatrix = systemMatrix;

      for (size_t i = 0; i < layerSolvers.size() - 1; i++)
      {
        almost::BuildInterpolationMatrix(layerSolvers[i].systemMatrix, layerSolvers[i].valuesCount, layerSolvers[i].interpolationMatrix, storage);
        layerSolvers[i].restrictionMatrix.BuildFromTransposed(layerSolvers[i].interpolationMatrix, storage);

        layerSolvers[i + 1].Resize(layerSolvers[i].restrictionMatrix.GetRowsCount());
        layerSolvers[i].tmpMatrix.BuildFromSparseProduct(layerSolvers[i].restrictionMatrix, layerSolvers[i].systemMatrix, storage);
        layerSolvers[i + 1].systemMatrix.BuildFromSparseProduct(layerSolvers[i].tmpMatrix, layerSolvers[i].interpolationMatrix, storage);
      }
    }

    template<typename StorageType>
    void Solve(const ValueType *rightSide, size_t iterationsCount, StorageType &storage)
    {
      std::copy_n(rightSide, layerSolvers[0].valuesCount, layerSolvers[0].rightSide.data());
      for (size_t i = 1; i < layerSolvers.size(); i++)
      {
        layerSolvers[i - 1].ResetValues();
        layerSolvers[i - 1].IterateGaussSeidel(iterationsCount / 2, storage);
        layerSolvers[i - 1].BuildResiduals(storage);
        BuildDenseVectorProduct(layerSolvers[i - 1].restrictionMatrix, layerSolvers[i - 1].residuals.data(), layerSolvers[i].rightSide.data(), storage);
      }
      layerSolvers.back().ResetValues();

      int p = 1;
      for (size_t i = 0; i < layerSolvers.size(); i++)
      {
        size_t layerIndex = layerSolvers.size() - 1 - i;
        layerSolvers[layerIndex].IterateGaussSeidel(iterationsCount / 2, storage);

        if (layerIndex > 0)
        {
          std::vector<ValueType> valuesDelta;
          valuesDelta.resize(layerSolvers[layerIndex - 1].valuesCount);

          almost::BuildDenseVectorProduct(layerSolvers[layerIndex - 1].interpolationMatrix, layerSolvers[layerIndex].values.data(), valuesDelta.data(), storage);
          for (size_t valueIndex = 0; valueIndex < layerSolvers[layerIndex - 1].valuesCount; valueIndex++)
            layerSolvers[layerIndex - 1].values[valueIndex] += valuesDelta[valueIndex];
        }
      }
    }

    ValueType *GetValues()
    {
      return layerSolvers[0].values.data();
    }

  private:
    std::vector<MultigridLayer<IndexType, ValueType>> layerSolvers;
  };
}