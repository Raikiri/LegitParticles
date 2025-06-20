#pragma once
#include <vector>

//#define INLINE __forceinline
#define INLINE

namespace almost
{
  template<typename DimTypeT>
  struct Index
  {
    using DimType = DimTypeT;
    using SelfType = Index<DimTypeT>;
    /*Index(const Index &other)
    {
      this->value = other.value;
      this->dim = other.dim;
    }*/
    Index(size_t _value) :
      value(_value)
    {
    }
    INLINE SelfType operator =(size_t value)
    {
      this->value = value;
    }
    INLINE SelfType &operator ++()
    {
      value++;
      return *this;
    }
    INLINE SelfType &operator *()
    {
      return *this;
    }
    INLINE bool operator != (const Index &other)
    {
      return value != other.value;
    }
    size_t value;
  };

  template<typename ...DimTypes>
  struct Space;

  template<class DimType>
  struct Space<DimType>
  {
    using SelfType = Space<DimType>;

    Space(DimType dim) : dim(dim)
    {
    }

    constexpr static size_t GetSizeStatic()
    {
      return DimType::size;
    }
    size_t GetSize()
    {
      return dim.size;
    }
    DimType dim;
  };

  template<class DimType, class ...DimTypes>
  struct Space<DimType, DimTypes...> : public Space<DimTypes...>
  {
    using SelfType = Space<DimType, DimTypes...>;
    using ParentType = Space<DimTypes...>;

    Space(DimType dim, DimTypes... dims) : dim(dim), ParentType(dims...)
    {
    }

    constexpr static size_t GetSizeStatic()
    {
      return DimType::size * ParentType::GetSizeStatic();
    }
    size_t GetSize()
    {
      return dim.size * this->ParentType::GetSize();
    }
    DimType dim;
  };

  template<typename IndexType>
  constexpr INLINE size_t GetStrideStatic(const IndexType &)
  {
    return IndexType::DimType::size;
  }

  template<typename IndexType, typename ...IndexTypes>
  constexpr INLINE size_t GetStrideStatic(const IndexType &, const IndexTypes &...)
  {
    return IndexType::DimType::size * GetStrideStatic<IndexTypes...>();
  }
    
  
  template<class IndexType>
  constexpr INLINE size_t GetLinearOffsetStatic(const IndexType& index)
  {
    return index.value;
  }
    
  template<class IndexType, class ...IndexTypes, typename = void>
  constexpr INLINE size_t GetLinearOffsetStatic(const IndexType &index, const IndexTypes& ...indices)
  {
    return index.value * GetStrideStatic<IndexTypes...>(indices...) + GetLinearOffsetStatic<IndexTypes...>(indices...);
  }
    
  template<typename IndexType>
  size_t GetLinearOffset(size_t &stride, IndexType index)
  {
    size_t resOffset = index.value;
    stride = index.dim.size;
    return resOffset;
  }

  template<typename IndexType, typename ...IndexTypes>
  size_t GetLinearOffset(size_t &stride, IndexType index, IndexTypes ...indices)
  {
    size_t resOffset = GetLinearOffset(stride, indices...) + index.value * stride;
    stride *= index.dim.size;
    return resOffset;
  }



  template<typename ValueTypeT, typename ...DimTypes>
  struct StaticTensor
  {
    using SpaceType = Space<DimTypes...>;
    using ValueType = ValueTypeT;

    template<typename ...IndexTypes>
    INLINE ValueType &Get(const Index<DimTypes> &...indices)
    {
      return data[GetLinearOffsetStatic<IndexTypes...>(indices...)];
    }
    template<typename ...IndexTypes>
    INLINE const ValueType& Get(const Index<DimTypes> &...indices) const
    {
      return data[GetLinearOffsetStatic<IndexTypes...>(indices...)];
    }

    ValueType data[SpaceType::GetSizeStatic()];
  };
  template<typename ValueTypeT, typename ...DimTypes>
  struct DynamicTensor
  {
    using SpaceType = Space<DimTypes...>;
    using ValueType = ValueTypeT;

    DynamicTensor(DimTypes... dimTypes)
    {
      data.resize(Space<DimTypes...>(dimTypes...).GetSize());
    }

    template<typename ...IndexTypes>
    INLINE ValueType &Get(Index<DimTypes> ...indices)
    {
      size_t stride;
      return data[GetLinearOffset<IndexTypes...>(stride, indices...)];
    }

    template<typename ...IndexTypes>
    INLINE const ValueType& Get(Index<DimTypes> ...indices) const
    {
      size_t stride;
      return data[GetLinearOffset<IndexTypes...>(stride, indices...)];
    }

    std::vector<ValueType> data;
  };

  template<typename ValueType, typename ...DimTypes>
  DynamicTensor<ValueType, DimTypes...> MakeDynamicTensor(DimTypes...dimTypes)
  {
    return DynamicTensor<ValueType, DimTypes...>(dimTypes...);
  }
  template<typename ValueType, typename ...DimTypes>
  constexpr StaticTensor<ValueType, DimTypes...> MakeStaticTensor()
  {
    return StaticTensor<ValueType, DimTypes...>();
  }


  template<typename DimensionName, size_t sizeT>
  struct StaticDimension
  {
    using SelfType = StaticDimension<DimensionName, sizeT>;
    using IndexType = Index<SelfType>;
    const static size_t size = sizeT;
    constexpr IndexType begin()
    {
      return IndexType(0);
    }
    constexpr IndexType end()
    {
      return IndexType(sizeT);
    }
  };

  template<typename DimensionName>
  struct DynamicDimension : public DimensionName
  {
    using SelfType = DynamicDimension<DimensionName>;
    using IndexType = Index<SelfType>;

    DynamicDimension(size_t size)
    {
      this->size = size;
    }
    IndexType begin()
    {
      return IndexType(*this, 0);
    }
    IndexType end()
    {
      return IndexType(*this, size);
    }
    size_t size;
  };
}
#undef INLINE

/*void UsageSample()
{
    struct Values {};
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
          valuesInBasis.Get(i, k) += valueToNullspace.Get(i, j) * nullspaceToBasis.Get(j, k);
}*/
