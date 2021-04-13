#pragma once
#include <algorithm>
#include <assert.h>
#include <vector>
namespace almost
{
  template<typename T>
  struct BasicSpace
  {
    using ValueType0 = T;
    using ValueType1 = T;
    using ResultType = T;
    using ScalarType = float;
    static T Mul(const ValueType0& val0, const T& val1)
    {
      return val0 * val1;
    }
    static ValueType0 Divide(const ResultType& val0, const T& val1)
    {
      return val0 / val1;
    }
    static ScalarType SpectralRadius(ValueType0 center, ValueType0 val0, ValueType0 val1)
    {
      return sqrt(center * center / abs(val0 * val1));
    }
  };
  template<typename RowDimensionT, typename ColumnDimensionT, typename ValueTypeT>
  struct SparseMatrix
  {
    using RowDimension = RowDimensionT;
    using ColumnDimension = ColumnDimensionT;
    using RowIndexType = typename RowDimension::IndexType;
    using ColumnIndexType = typename ColumnDimension::IndexType;
    using ValueType = ValueTypeT;
    using SelfType = SparseMatrix<RowDimension, ColumnDimension, ValueType>;
    struct RowTerm
    {
      ColumnIndexType columnIndex;
      ValueType value;
    };

    /*void AddRow(RowIndexType rowIndex, const ColumnIndexType *columnIndices, const ValueType *values, size_t termsCount)
    {
      Row newRow;
      newRow.rowIndex = rowIndex;
      newRow.termsStart = rowTerms.size();
      newRow.termsCount = termsCount;

      rowTerms.resize(rowTerms.size() + termsCount);
      for (size_t termIndex = 0; termIndex < termsCount; termIndex++)
      {
        rowTerms[newRow.termsStart + termIndex].columnIndex = columnIndices[termIndex];
        rowTerms[newRow.termsStart + termIndex].value = values[termIndex];
      }
      rows.push_back(newRow);
    }*/

    void AppendRow(RowIndexType rowIndex)
    {
      Row newRow;
      newRow.rowIndex = rowIndex;
      newRow.termsStart = rowTerms.size();
      newRow.termsCount = 0;
      rows.push_back(newRow);
    }


    void AppendTerm(ColumnIndexType columnIndex, ValueType value)
    {
      auto& lastRow = rows.back();
      lastRow.termsCount++;

      RowTerm newTerm;
      newTerm.columnIndex = columnIndex;
      newTerm.value = value;
      rowTerms.push_back(newTerm);
    }

    void AppendTerms(RowTerm* newRowTerms, size_t termsCount)
    {
      auto& lastRow = rows.back();
      lastRow.termsCount += termsCount;      
      std::copy_n(newRowTerms, termsCount, std::back_inserter(rowTerms));
    }

    size_t GetRowsCount() const
    {
      return rows.size();
    }

    RowIndexType GetRowIndex(size_t rowNumber) const
    {
      return rows[rowNumber].rowIndex;
    }
    size_t GetRowTermsCount(size_t rowNumber) const
    {
      return rows[rowNumber].termsCount;
    }
    RowTerm* GetRowTerms(size_t rowNumber)
    {
      size_t termsStart = rows[rowNumber].termsStart;
      return rowTerms.data() + termsStart;
    }
    const RowTerm* GetRowTerms(size_t rowNumber) const
    {
      size_t termsStart = rows[rowNumber].termsStart;
      return rowTerms.data() + termsStart;
    }

    bool CheckSortedIndices() const
    {
      //return true;
      bool isOk = true;
      for (size_t rowNumber = 0; rowNumber < rows.size(); rowNumber++)
      {
        if (rowNumber < rows.size() - 1)
          isOk = isOk && (rows[rowNumber].rowIndex < rows[rowNumber + 1].rowIndex);

        auto& row = rows[rowNumber];
        for (size_t columnNumber = 0; columnNumber < row.termsCount - 1; columnNumber++)
        {
          isOk = isOk && (rowTerms[row.termsStart + columnNumber].columnIndex < rowTerms[row.termsStart + columnNumber + 1].columnIndex);
        }
      }
      assert(isOk);
      return isOk;
    }

    /*void SortIndices()
    {
      std::sort(rows.begin(), rows.end(),
        [](const Row& left, const Row& right) -> bool
        {
          assert(left.rowIndex < right.rowIndex || right.rowIndex < left.rowIndex); // can't have 2 identical row indices
          return left.rowIndex < right.rowIndex;
        });

      for (auto row : rows)
      {
        std::sort(rowTerms.begin() + row.termsStart, rowTerms.begin() + row.termsStart + row.termsCount,
          [](const RowTerm& left, const RowTerm& right) -> bool
          {
            assert(left.columnIndex < right.columnIndex || right.columnIndex < left.columnIndex); // can't have 2 identical row indices
            return left.columnIndex < right.columnIndex;
          });
      }
    };
    SelfType &BuildSorted(const SelfType &src)
    {
      this->rowDimension = src.rowDimension;
      this->columnDimension = src.columnDimension;
      this->rows = src.rows;
      this->rowTerms.resize(src.rowTerms.size());

      size_t rowTermsOffset = 0;
      for (auto &row : rows)
      {
        std::copy(src.rowTerms.begin() + row.termsStart, src.rowTerms.begin() + row.termsStart + row.termsCount, this->rowTerms.begin() + rowTermsOffset);
        row.termsStart = rowTermsOffset;
        rowTermsOffset += row.termsCount;
      }
      return *this;
    }*/

    struct Element
    {
      RowIndexType rowIndex;
      ColumnIndexType columnIndex;
      ValueType value;
    };

    static void SortElements(Element *elements, size_t count)
    {
      std::sort(elements, elements + count,
        [](const Element& left, const Element& right) -> bool
        {
          if (left.rowIndex < right.rowIndex) return true;
          if (right.rowIndex < left.rowIndex) return false;
          return left.columnIndex < right.columnIndex;
        });
    }
    SelfType &BuildFromSortedElements(const Element *elements, size_t count, RowDimension rowDimension, ColumnDimension columnDimension)
    {
      BuildEmpty(rowDimension, columnDimension);
      for (size_t index = 0; index < count; index++)
      {
        /*if (elements[index].value == ValueType(0))
          continue;*/

        bool sameRow = (rows.size() > 0) && RowsEqual(rows.back().rowIndex, elements[index].rowIndex);
        bool sameColumn = (rowTerms.size() > 0) && ColumnsEqual(rowTerms.back().columnIndex, elements[index].columnIndex);

        if (sameRow && sameColumn)
        {
          rowTerms.back().value += elements[index].value;
        }
        else
        {
          if (sameRow)
          {
            rows.back().termsCount++;
          }else
          {
            Row newRow;
            newRow.rowIndex = elements[index].rowIndex;
            newRow.termsCount = 1;
            newRow.termsStart = rowTerms.size();
            rows.push_back(newRow);
          }

          RowTerm newTerm;
          newTerm.columnIndex = elements[index].columnIndex;
          newTerm.value = elements[index].value;
          rowTerms.push_back(newTerm);
        }
      }
      assert(CheckSortedIndices());
      return *this;
    }

    template<typename StorageType>
    SelfType &BuildFromTransposed(const SparseMatrix<ColumnDimension, RowDimension, ValueType> &srcMatrix, StorageType &storage)
    {
      auto elements = storage.template GetHandle<std::vector<Element> >( );
      elements.Get().clear();
      for (const auto &row : srcMatrix.rows)
      {
        for (size_t termIndex = 0; termIndex < row.termsCount; termIndex++)
        {
          Element newElement;
          newElement.columnIndex = row.rowIndex;
          const auto &term = srcMatrix.rowTerms[row.termsStart + termIndex];
          newElement.rowIndex = term.columnIndex;
          newElement.value = term.value;
          elements.Get().push_back(newElement);
        }
      }

      SelfType::SortElements(elements.Get().data(), elements.Get().size());
      return BuildFromSortedElements(elements.Get().data(), elements.Get().size(), srcMatrix.columnDimension, srcMatrix.rowDimension);
    }

    template<typename Space, typename ValueType0, typename ValueType1, typename CommonDimension, typename StorageType>
    SelfType &BuildFromDenseProduct(const SparseMatrix<RowDimension, CommonDimension, ValueType0>& matrix0, const SparseMatrix<ColumnDimension, CommonDimension, ValueType1> &matrix1Transposed, StorageType &storage)
    {
      assert(matrix0.columnDimension.size == matrix1Transposed.columnDimension.size);

      BuildEmpty(matrix0.rowDimension, matrix1Transposed.rowDimension);
      auto elements = storage.template GetHandle<std::vector<Element>>();
      elements.Get().clear();

      size_t rowNumber0 = 0;
      size_t rowNumber1 = 0;
      for (size_t rowNumber0 = 0; rowNumber0 < matrix0.rows.size(); rowNumber0++)
      {
        for (size_t rowNumber1 = 0; rowNumber1 < matrix1Transposed.rows.size(); rowNumber1++)
        {
          auto &row0 = matrix0.rows[rowNumber0];
          auto &row1 = matrix1Transposed.rows[rowNumber1];

          Element newElement;
          newElement.rowIndex = row0.rowIndex;
          newElement.columnIndex = row1.rowIndex;

          newElement.value = ValueType(0.0f);

          size_t termNumber0 = 0;
          size_t termNumber1 = 0;
          while (termNumber0 < row0.termsCount && termNumber1 < row1.termsCount)
          {
            auto &term0 = matrix0.rowTerms[row0.termsStart + termNumber0];
            auto &term1 = matrix1Transposed.rowTerms[row1.termsStart + termNumber1];
            if (term0.columnIndex < term1.columnIndex)
              termNumber0++;
            else
            if (term1.columnIndex < term0.columnIndex)
              termNumber1++;
            else
            {
              newElement.value += Space::Get(term0.value, term1.value);
              termNumber0++;
              termNumber1++;
            }
          }
          if(newElement.value != ValueType(0))
            elements.Get().push_back(newElement);
        }
      }

      return BuildFromSortedElements(elements.Get().data(), elements.Get().size(), matrix0.rowDimension, matrix1Transposed.rowDimension);
    }
    
    template<typename StorageType>
    SelfType &BuildFromSum(const SelfType& matrix0, float mult0, const SelfType& matrix1, float mult1, StorageType &storage)
    {
      assert(matrix0.columnDimension.size == matrix1.columnDimension.size);
      assert(matrix0.rowDimension.size == matrix1.rowDimension.size);

      BuildEmpty(matrix0.rowDimension, matrix0.columnDimension);
      auto elementsHandle = storage.template GetHandle < std::vector < Element >> ();
      auto elements = elementsHandle.Get();
      elements.clear();

      
      for (size_t rowNumber = 0; rowNumber < matrix0.rows.size(); rowNumber++)
      {
        RowIndexType rowIndex = matrix0.GetRowIndex(rowNumber);
        size_t termsCount = matrix0.GetRowTermsCount(rowNumber);
        const RowTerm* terms = matrix0.GetRowTerms(rowNumber);
        for (size_t termNumber = 0; termNumber < termsCount; termNumber++)
        {
          auto& term = terms[termNumber];
          Element newElement;
          newElement.rowIndex = rowIndex;
          newElement.columnIndex = term.columnIndex;
          newElement.value = term.value * mult0;
          elements.push_back(newElement);
        }
      }

      for (size_t rowNumber = 0; rowNumber < matrix1.rows.size(); rowNumber++)
      {
        RowIndexType rowIndex = matrix1.GetRowIndex(rowNumber);
        size_t termsCount = matrix1.GetRowTermsCount(rowNumber);
        const RowTerm* terms = matrix1.GetRowTerms(rowNumber);
        for (size_t termNumber = 0; termNumber < termsCount; termNumber++)
        {
          auto& term = terms[termNumber];
          Element newElement;
          newElement.rowIndex = rowIndex;
          newElement.columnIndex = term.columnIndex;
          newElement.value = term.value * mult1;
          elements.push_back(newElement);
        }
      }

      SortElements(elements.data(), elements.size());
      return BuildFromSortedElements(elements.data(), elements.size(), matrix0.rowDimension, matrix0.columnDimension);

      /*auto columnIndicesHandle0 = storage.template GetHandle<std::vector<ColumnIndexType>>();
      auto columnTermsHandle0 = storage.template GetHandle<std::vector<ValueType>>();
      auto columnIndicesHandle1 = storage.template GetHandle<std::vector<ColumnIndexType>>();
      auto columnTermsHandle1 = storage.template GetHandle<std::vector<ValueType>>();
      auto& columnIndices0 = columnIndicesHandle0.Get();
      auto& columnTerms0 = columnTermsHandle0.Get();
      auto& columnIndices1 = columnIndicesHandle1.Get();
      auto& columnTerms1 = columnTermsHandle1.Get();
      indices.clear();
      values.clear();


      size_t rowNumber0 = 0;
      size_t rowNumber1 = 0;
      while (rowNumber0 < matrix0.rows.size() && rowNumber1 < matrix1.rows.size();)
      {
        auto& row0 = matrix0.rows[rowNumber0];
        auto& row1 = matrix1.rows[rowNumber1];

        if (row0.rowIndex < row1.rowIndex)
        {
          size_t termsCount = matrix0.GetRowTermsCount(rowNumber0);
          columnIndices0.resize(termsCount);
          columnTerms0.resize(termsCount);

          matrix0.GetRowTerms(rowNumber0, columnIndices0(), columnTerms0.data());
          AddRow(row0.rowIndex, columnIndices0.data(), columnTerms0.data(), termsCount);
          rowNumber0++;
        }
        else
        if (row1.rowIndex < row0.rowIndex)
        {
          size_t termsCount = matrix1.GetRowTermsCount(rowNumber1);
          columnIndices1.resize(termsCount);
          columnTerms1.resize(termsCount);

          matrix1.GetRowTerms(rowNumber1, columnIndices1.data(), columnTerms1.data());
          AddRow(row1.rowIndex, columnIndices1.data(), columnTerms1.data(), termsCount);
          rowNumber1++;
        }
        else
        {
          newElement.value += Space::Get(term0.value, term1.value);
          termNumber0++;
          termNumber1++;
        }
      }

      return BuildFromSortedElements(elements.Get().data(), elements.Get().size(), matrix0.rowDimension, matrix1Transposed.rowDimension);*/
    }


    template<typename Space, typename CommonDimension, typename ValueType0, typename ValueType1, typename StorageType>
    SelfType &BuildFromSparseProduct(const SparseMatrix<RowDimension, CommonDimension, ValueType0>& matrix0, const SparseMatrix<CommonDimension, ColumnDimension, ValueType1>& matrix1, StorageType& storage)
    {
      assert(matrix0.columnDimension.size == matrix1.rowDimension.size);
      BuildEmpty(matrix0.rowDimension, matrix1.columnDimension);
      auto elementsHandle = storage.template GetHandle<std::vector<Element>>();
      auto& elements = elementsHandle.Get();
      elements.clear();
      //matrix0.CheckSortedIndices();
      //matrix1.CheckSortedIndices();
      for (size_t rowNumber = 0; rowNumber < matrix0.rows.size(); rowNumber++)
      {
        auto &row = matrix0.rows[rowNumber];
        size_t rowStart = elements.size();
        for (size_t termNumber = 0; termNumber < row.termsCount; termNumber++)
        {
          auto &term = matrix0.rowTerms[row.termsStart + termNumber];
          size_t otherRowNumber = matrix1.FindRowNumber(term.columnIndex);
          if (otherRowNumber == size_t(-1))
            continue;
          auto &otherRow = matrix1.rows[otherRowNumber];


          for (size_t otherTermNumber = 0; otherTermNumber < otherRow.termsCount; otherTermNumber++)
          {
            auto &otherTerm = matrix1.rowTerms[otherRow.termsStart + otherTermNumber];
            auto &otherColumnIndex = otherTerm.columnIndex;
            Element newElement;
            newElement.rowIndex = row.rowIndex;
            newElement.columnIndex = otherColumnIndex;
            newElement.value = Space::Mul(term.value, otherTerm.value);
            //if (newElement.value != ProductType(0))
            elements.push_back(newElement);
          }
        }
        SortElements(elements.data() + rowStart, elements.size() - rowStart);
      }
      //SortElements(elements.Get().data(), elements.Get().size());
      return BuildFromSortedElements(elements.data(), elements.size(), matrix0.rowDimension, matrix1.columnDimension);
    }

    SelfType& BuildFromDiag(ValueType diag, RowDimension rowDimension)
    {
      BuildEmpty(rowDimension, rowDimension);
      for (RowIndexType rowIndex = 0; rowIndex < rowDimension.size; rowIndex++)
      {
        AppendRow(rowIndex);
        AppendTerm(rowIndex, diag);
      }
      return *this;
    }
    size_t FindRowNumber(RowIndexType rowIndex) const
    {
      Row tmpRow;
      tmpRow.rowIndex = rowIndex;
      auto it = std::lower_bound(rows.begin(), rows.end(), tmpRow, [](const Row &left, const Row &right) -> bool
      {
        return left.rowIndex < right.rowIndex;
      });
      if (it != rows.end() && RowsEqual(it->rowIndex, rowIndex))
        return std::distance(rows.begin(), it);
      else
        return size_t(-1);
    }
    bool FindTerm(size_t rowNumber, ColumnIndexType columnIndex, ValueType &value) const
    {
      const auto &row = rows[rowNumber];

      RowTerm tmpTerm;
      tmpTerm.columnIndex = columnIndex;
      auto termsStart = rowTerms.begin() + row.termsStart;
      auto termsEnd = termsStart + row.termsCount;
      auto it = std::lower_bound(termsStart, termsEnd, tmpTerm, [](const RowTerm &left, const RowTerm &right) -> bool
      {
        return left.columnIndex < right.columnIndex;
      });
      if (it != termsEnd && ColumnsEqual(it->columnIndex, columnIndex))
      {
        value = it->value;
        return true;
      }
      return false;
    }

    void BuildEmpty(RowDimension rowDimension, ColumnDimension columnDimension)
    {
      this->rowDimension = rowDimension;
      this->columnDimension = columnDimension;
      rows.clear();
      rowTerms.clear();
    }
    static bool RowsEqual(RowIndexType left, RowIndexType right)
    {
      return !(left < right) && !(right < left);
    };
    static bool ColumnsEqual(ColumnIndexType left, ColumnIndexType right)
    {
      return !(left < right) && !(right < left);
    };
    struct Row
    {
      RowIndexType rowIndex;
      size_t termsStart;
      size_t termsCount;
    };
    std::vector<Row> rows;

    std::vector<RowTerm> rowTerms;
    RowDimension rowDimension;
    ColumnDimension columnDimension;
  };
}