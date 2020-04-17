#define INLINE __forceinline


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
    Index(DimType _dim, size_t _value) :
      dim(_dim), value(_value)
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
    INLINE bool operator != (Index &other)
    {
      return value != other.value;
    }
    size_t value;
    DimType dim;
  };

  template<typename DimType>
  Index<DimType> MakeIndex(DimType dim, size_t value = 0)
  {
    return Index<DimType>(dim, value);
  }

  template<class ...>
  struct CombinedIndex;

  template<class IndexType>
  struct CombinedIndex<IndexType>
  {
    using SelfType = CombinedIndex<IndexType>;
    CombinedIndex(IndexType &index) : index(index)
    {
      linearOffset = index.value;
      stride = index.dim.size;
      carry = 0;
    }
    INLINE SelfType &operator++()
    {
      linearOffset++;
      index.value++;
      carry = (index.value >= index.dim.size);
      index.value = carry ? 0 : index.value;
      return *this;
    }
    INLINE size_t GetLinearOffset() const
    {
      return linearOffset;
    }

    INLINE SelfType begin()
    {
      Init();
      return *this;
    }
    INLINE SelfType end()
    {
      return *this;
    }
    INLINE SelfType &operator *()
    {
      return *this;
    }
    INLINE bool operator != (const SelfType &other)
    {
      return carry == 0; //just used for iterating, not actual comparison
    }
  protected:
    void Init()
    {
      index.value = 0;
    }
    IndexType &index;
    size_t carry;
    size_t linearOffset;
    size_t stride;
  };

  /*template<class ...IndexTypes>
  CombinedIndex<IndexTypes...> CreateStaticIndex();*/

  template<class IndexType, class ...IndexTypes>
  struct CombinedIndex<IndexType, IndexTypes...> : public CombinedIndex<IndexTypes...>
  {
    using ParentType = CombinedIndex<IndexTypes...>;
    using SelfType = CombinedIndex<IndexType, IndexTypes...>;

    CombinedIndex(IndexType &index, IndexTypes&... indices) : index(index), ParentType(indices...)
    {
      this->linearOffset += this->ParentType::stride * index.value;
      this->stride *= index.dim.size; 
      this->carry = 0;
    }
    INLINE SelfType &operator++()
    {
      ParentType::operator++();
      index.value += this->ParentType::carry;
      carry = (index.value >= index.dim.size);
      index.value = carry ? 0 : index.value;
      return *this;
    }
    constexpr size_t GetSizeStatic()
    {
      return IndexType::DimType::size * ParentType::GetSizeStatic();;
    }
    size_t GetSize()
    {
      return index.dim.size * this->ParentType::GetSize();
    }
    INLINE SelfType &begin()
    {
      Init();
      return *this;
    }
    INLINE SelfType &end()
    {
      return *this;
    }
    INLINE SelfType &operator *()
    {
      return *this;
    }
    INLINE bool operator != (const SelfType &other)
    {
      return carry == 0; //just used for iterating, not actual comparison
    }
  protected:
    void Init()
    {
      index.value = 0;
      ParentType::Init();
    }
    IndexType &index;
    size_t carry;
  };
  
  template<typename ...IndexTypes>
  CombinedIndex<IndexTypes...> MakeCombinedIndex(IndexTypes&...indices)
  {
    return CombinedIndex<IndexTypes...>(indices...);
  }

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
  constexpr size_t GetStrideStatic(IndexType index)
  {
    return IndexType::DimType::size;
  }

  template<typename IndexType, typename ...IndexTypes>
  constexpr size_t GetStrideStatic(IndexType index, IndexTypes...indices)
  {
    return IndexType::DimType::size * GetStrideStatic<IndexTypes...>(indices);
  }
    

  template<class IndexType>
  constexpr size_t GetLinearOffsetStatic(IndexType index)
  {
    return index.value;
  }
    
  template<class IndexType, class ...IndexTypes>
  constexpr size_t GetLinearOffsetStatic(IndexType index, IndexTypes ...indices)
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
    using CombinedIndexType = CombinedIndex<Index<DimTypes>...>;
    using SpaceType = Space<DimTypes...>;
    using ValueType = ValueTypeT;

    template<typename ...IndexTypes>
    INLINE ValueType &Get(Index<DimTypes> ...indices)
    {
      return data[GetLinearOffsetStatic<IndexTypes...>(indices...)];
    }

    ValueType data[SpaceType::GetSizeStatic()];
  };
  template<typename ValueTypeT, typename ...DimTypes>
  struct DynamicTensor
  {
    using CombinedIndexType = CombinedIndex<Index<DimTypes>...>;
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

    std::vector<ValueType> data;
  };

  template<typename ValueType, typename ...DimTypes>
  DynamicTensor<ValueType, DimTypes...> MakeDynamicTensor(DimTypes...dimTypes)
  {
    return DynamicTensor<ValueType, DimTypes...>(dimTypes...);
  }
  template<typename ValueType, typename ...DimTypes>
  StaticTensor<ValueType, DimTypes...> MakeStaticTensor()
  {
    return StaticTensor<ValueType, DimTypes...>();
  }


  template<typename DimensionName, size_t sizeT>
  struct StaticDimension : public DimensionName
  {
    using SelfType = StaticDimension<DimensionName, sizeT>;
    using IndexType = Index<SelfType>;
    const static size_t size = sizeT;
    IndexType begin()
    {
      return IndexType(*this, 0);
    }
    IndexType end()
    {
      return IndexType(*this, sizeT);
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