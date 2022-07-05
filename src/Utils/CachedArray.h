#pragma once
#include <vector>
#include <utility>
#include <assert.h>
namespace almost
{

struct Reindex
{
  Reindex() {}
  Reindex(size_t srcIndex, size_t dstIndex)
  {
    this->srcIndex = srcIndex;
    this->dstIndex = dstIndex;
  }
  size_t srcIndex;
  size_t dstIndex;
};

template<typename TypeName>
struct Id
{
public:
  Id() {}
  Id(size_t _id) : id(_id) {}
  size_t id;
};


template<typename T, typename TypeName = T>
struct CachedArray
{
public:
  using Id = Id<TypeName>;
  T &GetById(Id elementId)
  {
    return elements[elementId.id];
  }
  T &GetByIndex(size_t elementIndex)
  {
    return elements[aliveIds[elementIndex].id];
  }
  T &operator[](size_t elementIndex)
  {
    return GetByIndex(elementIndex);
  }
  T &operator[](Id elementId)
  {
    return GetById(elementId);
  }
  Id GetId(size_t elementIndex)
  {
    return aliveIds[elementIndex];
  }
  size_t GetElementsCount()
  {
    return aliveIds.size();
  }
  size_t GetAllocatedCount()
  {
    return elements.size();
  }
  Id Add(bool &isCreated)
  {
    if (freeIds.size() > 0)
    {
      Id freeElementId = freeIds.back();
      freeIds.pop_back();
      aliveIds.push_back(freeElementId);
      states[freeElementId.id] = State::Alive;
      isCreated = false;
      return freeElementId;
    }
    else
    {
      aliveIds.push_back(elements.size());
      elements.push_back(T());
      states.push_back(State::Alive);
      isCreated = true;
      return Id(elements.size() - 1);
    }
  }
  Id Add()
  {
    bool isCreated;
    return Add(isCreated);
  }
  Id Add(const T &newbie)
  {
    bool isCreated;
    return Add(newbie, isCreated);
  }
  Id Add(const T &newbie, bool &isCreated)
  {
    if(freeIds.size() > 0)
    {
      Id freeElementId = freeIds.back();
      freeIds.pop_back();
      aliveIds.push_back(freeElementId);
      states[freeElementId.id] = State::Alive;
      elements[freeElementId.id] = newbie;
      isCreated = false;
      return freeElementId;
    }else
    {
      aliveIds.push_back(elements.size());
      elements.push_back(newbie);
      states.push_back(State::Alive);
      isCreated = true;
      return elements.size() - 1;
    }
  }
  void RemoveById(Id elementId)
  {
    if (states[elementId.id] == State::Alive)
      states[elementId.id] = State::Releasing;
  }
  void RemoveByIndex(size_t elementIndex)
  {
    if(states[aliveIds[elementIndex].id] == State::Alive)
      states[aliveIds[elementIndex].id] = State::Releasing;
  }
  bool IsIdAlive(Id elementId)
  {
    if (elementId.id < states.size())
    {
      return states[elementId.id] == State::Alive;
    }
    else
      return false;
  }
  void Update()
  {
    aliveIds.clear();
    for(size_t elementId = 0; elementId < elements.size(); elementId++)
    {
      if(states[elementId] == State::Releasing)
      {
        states[elementId] = State::Dead;
        freeIds.push_back(Id(elementId));
      }
      if(states[elementId] == State::Alive)
      {
        aliveIds.push_back(Id(elementId));
      }
    }
  }
private:
  struct State
  {
    enum Types
    {
      Alive,
      Releasing,
      Dead
    };
  };
  std::vector<Id> freeIds;
  std::vector<Id> aliveIds;
  std::vector<typename State::Types> states;
  std::vector<T> elements;
};

template<typename T, typename ContainerType>
struct LinearIterator
{
  using SelfType = LinearIterator<T, ContainerType>;
  LinearIterator(ContainerType* container, size_t index) : container(container), index(index) {}
  T& operator *()
  {
    return container->GetByIndex(index);
  }
  bool operator != (const SelfType& other)
  {
    return this->index != other.index; //just used for iterating, not actual comparison
  }
  SelfType& operator ++()
  {
    index++;
    return *this;
  }
  ContainerType* container;
  size_t index;
};

template<typename T, typename TypeName = T>
struct LinearCachedArray
{
public:
  using SelfType = LinearCachedArray<T, TypeName>;
  using Id = Id<TypeName>;
  using Iterator = LinearIterator<T, SelfType>;
  LinearCachedArray()
  {
    this->elementsAdded   = false;
    this->elementsRemoved = false;
  }
  T &GetById(Id elementId)
  {
    return elements[idToIndex[elementId.id]];
  }
  T &GetByIndex(size_t elementIndex)
  {
    return elements[elementIndex];
  }
  T *GetElements()
  {
    return elements.data();
  }
  size_t GetElementsCount()
  {
    return elements.size();
  }
  Id GetId(size_t elementIndex)
  {
    return indexToId[elementIndex];
  }
  size_t GetIndex(Id elementId)
  {
    return idToIndex[elementId.id];
  }

  size_t GetAllocatedCount()
  {
    return elements.capacity();
  }

  Id Add(bool &isCreated)
  {
    return Add(T(), isCreated);
  }
  Id Add()
  {
    bool isCreated;
    return Add(T(), isCreated);
  }
  Id Add(T &&newbie)
  {
    bool isCreated;
    return Add(std::forward<T>(newbie), isCreated);
  }
  Id Add(T &&newbie, bool &isCreated)
  {
    this->elementsAdded = true;
    Id newId = Id(-1);
    if (freeIds.size() > 0)
    {
      newId = freeIds.back();
      freeIds.pop_back();

      idToIndex[newId.id] = elements.size();

      isCreated = false;
    }
    else
    {
      newId = elements.size();
      size_t newIndex = elements.size();

      assert(idToIndex.size() == elements.size());
      idToIndex.push_back(newIndex);

      isCreated = true;
    }
    assert(indexToId.size() == elements.size());
    indexToId.push_back(newId);

    assert(states.size() == elements.size());
    states.push_back(State::Alive);

    elements.emplace_back(std::forward<T>(newbie));
    return newId;
  }
  void RemoveById(Id elementId)
  {
    this->elementsRemoved = true;
    if (states[idToIndex[elementId.id]] == State::Alive)
      states[idToIndex[elementId.id]] = State::Releasing;
  }
  void RemoveByIndex(size_t elementIndex)
  {
    this->elementsRemoved = true;
    if (states[elementIndex] == State::Alive)
      states[elementIndex] = State::Releasing;
  }
  bool IsIdAlive(Id elementId)
  {
    if (elementId.id < idToIndex.size())
    {
      return states[idToIndex[elementId.id]] == State::Alive;
    }
    else
      return false;
  }
  bool IsModified()
  {
    return elementsAdded || elementsRemoved;
  }
  bool ElementsAdded()
  {
    return elementsAdded;
  }
  bool ElementsRemoved()
  {
    return elementsRemoved;
  }



  size_t GetReindicesCount()
  {
    return reindices.size();
  }
  Reindex *GetReindices()
  {
    return reindices.data();
  }

  void Update()
  {
    reindices.clear();
    if (IsModified())
    {
      for (size_t elementIndex = 0; elementIndex < elements.size();)
      {
        if (states[elementIndex] == State::Releasing)
        {
          size_t lastIndex = elements.size() - 1;
          Id lastId = indexToId[lastIndex];

          Id elementId = indexToId[elementIndex];

          freeIds.push_back(elementId);

          indexToId[lastIndex] = Id(-1);
          idToIndex[lastId.id] = elementIndex;

          indexToId[elementIndex] = Id(lastId);
          idToIndex[elementId.id] = -1;

          states[elementIndex] = states[lastIndex];
          states[lastIndex] = State::Dead;

          //elements[elementIndex].~T();
          elements[elementIndex] = std::move(elements[lastIndex]);
          reindices.push_back(Reindex(lastIndex, elementIndex));
          reindices.push_back(Reindex(elementIndex, size_t(-1)));

          indexToId.pop_back();
          states.pop_back();
          elements.pop_back();
        }
        else
        {
          elementIndex++;
        }
      }
    }
    elementsAdded = false;
    elementsRemoved = false;
  }

  Iterator begin()
  {
    return Iterator(this, 0);
  }
  Iterator end()
  {
    return Iterator(this, GetElementsCount());
  }
private:
  struct State
  {
    enum Types
    {
      Alive,
      Releasing,
      Dead
    };
  };
  std::vector<Id> freeIds;
  std::vector<size_t> idToIndex;
  std::vector<Id> indexToId;

  bool elementsAdded;
  bool elementsRemoved;


  std::vector<Reindex> reindices;

  std::vector<typename State::Types> states;
  std::vector<T> elements;
};

}