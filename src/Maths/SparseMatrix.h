#pragma once
#include <algorithm>
#include <assert.h>
#include <vector>
namespace almost
{
  template<typename RowIndexTypeT, typename ColumnIndexTypeT, typename ValueTypeT>
  struct SparseMatrix
  {
    using RowIndexType = RowIndexTypeT;
    using ColumnIndexType = ColumnIndexTypeT;
    using ValueType = ValueTypeT;
    using SelfType = SparseMatrix<RowIndexType, ColumnIndexType, ValueType>;
    void AddRow(RowIndexType rowIndex, const ColumnIndexType *columnIndices, const ValueType *values, size_t termsCount)
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
    void GetRowTerms(size_t rowNumber, ColumnIndexType *columnIndices, ValueType *values) const
    {
      for (size_t termIndex = 0; termIndex < rows[rowNumber].termsCount; termIndex++)
      {
        size_t termsStart = rows[rowNumber].termsStart;
        columnIndices[termIndex] = rowTerms[termsStart + termIndex].columnIndex;
        values[termIndex] = rowTerms[termsStart + termIndex].value;
      }
    }

    bool CheckSortedIndices()
    {
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

    void SortIndices()
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
    SelfType &BuildSorted(const SparseMatrix<RowIndexType, ColumnIndexType, ValueType> &src)
    {
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
    }

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
    SelfType &BuildFromSortedElements(const Element *elements, size_t count)
    {
      rows.clear();
      rowTerms.clear();
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
      return *this;
    }

    template<typename StorageType>
    SelfType &BuildFromTransposed(const SparseMatrix<ColumnIndexType, RowIndexType, ValueType> &srcMatrix, StorageType &storage)
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
      return BuildFromSortedElements(elements.Get().data(), elements.Get().size());
    }

    template<typename CommonIndexType, typename StorageType>
    SelfType &BuildFromDenseProduct(const SparseMatrix<RowIndexType, CommonIndexType, ValueType>& matrix0, const SparseMatrix<ColumnIndexType, CommonIndexType, ValueType> &matrix1Transposed, StorageType &storage)
    {
      Clear();
      auto elements = storage.template GetHandle<std::vector<Element>>();
      elements.Get().clear();

      size_t rowNumber0 = 0;
      size_t rowNumber1 = 0;
      for (size_t rowNumber1 = 0; rowNumber1 < matrix1Transposed.rows.size(); rowNumber1++)
      {
        for(size_t rowNumber0 = 0; rowNumber0 < matrix0.rows.size(); rowNumber0++)
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
              newElement.value += term0.value * term1.value;
              termNumber0++;
              termNumber1++;
            }
          }
          if(newElement.value != ValueType(0))
            elements.Get().push_back(newElement);
        }
      }

      return BuildFromSortedElements(elements.Get().data(), elements.Get().size());
    }

    template<typename CommonIndexType, typename ValueType0, typename ValueType1, typename StorageType>
    SelfType &BuildFromSparseProduct(const SparseMatrix<RowIndexType, CommonIndexType, ValueType0>& matrix0, const SparseMatrix<CommonIndexType, ColumnIndexType, ValueType1>& matrix1, StorageType& storage)
    {
      Clear();
      auto elements = storage.template GetHandle<std::vector<Element>>();
      elements.Get().clear();

      for (size_t rowNumber = 0; rowNumber < matrix0.rows.size(); rowNumber++)
      {
        auto &row = matrix0.rows[rowNumber];
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
            newElement.value = term.value * otherTerm.value;
            //if (newElement.value != ProductType(0))
            elements.Get().push_back(newElement);
          }
        }
      }

      return BuildFromSortedElements(elements.Get().data(), elements.Get().size());
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

    void Clear()
    {
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

    struct RowTerm
    {
      ColumnIndexType columnIndex;
      ValueType value;
    };
    std::vector<RowTerm> rowTerms;
  };
}