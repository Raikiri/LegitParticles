#pragma once
#include <vector>
namespace almost
{
  template<typename T> struct Overload {};

  template<typename ...Rest>
  struct StackStorageData;

  template<typename First, typename ...Rest>
  struct StackStorageData<First, Rest...> : public StackStorageData<Rest...>
  {
    using StackStorageData<Rest...>::Push;
    using StackStorageData<Rest...>::Pop;
    First& Push(Overload<First>)
    {
      if (values.size() < usedCount + 1)
      {
        values.resize(usedCount + 1);
        values[usedCount] = std::make_unique<First>();
      }
      return *values[usedCount++];
    }
    void Pop(Overload<First>)
    {
      usedCount--;
    }
    std::vector<std::unique_ptr<First>> values;
    size_t usedCount = 0;
  };

  template<>
  struct StackStorageData<>
  {
    void Push() {}
    void Pop() {}
  };

  template<typename T, typename StorageType>
  struct StackHandle
  {
    StackHandle(T& value_, StorageType& stackStorage) :
      stackStorage(stackStorage),
      value(value_) {}
    ~StackHandle() { stackStorage.Pop(Overload<T>()); }
    T& Get() { return value; }
    T& value;
    StorageType& stackStorage;
  };

  template<typename ...Types>
  struct StackStorage : public StackStorageData<Types...>
  {
    using StackStorageData<Types...>::Push;
    using StackStorageData<Types...>::Pop;

    using SelfType = StackStorage<Types...>;
    template<typename T>
    StackHandle<T, SelfType> GetHandle()
    {
      return StackHandle<T, SelfType>(Push(Overload<T>()), *this);
    }
  };
}