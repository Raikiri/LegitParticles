#pragma once

#include "json/json.h"
#include "../../Utils/GlmInclude.h"

#include <iostream>
#include <algorithm>
namespace almost
{
  struct PackingResult
  {
    std::vector<glm::ivec2> rectPositions;
    glm::ivec2 size;
  };

  bool TestRects(glm::ivec2 pos0, glm::ivec2 size0, glm::ivec2 pos1, glm::ivec2 size1, int padding)
  {
    return
      (pos1.x > pos0.x + size0.x + padding) ||
      (pos0.x > pos1.x + size1.x + padding) ||
      (pos1.y > pos0.y + size0.y + padding) ||
      (pos0.y > pos1.y + size1.y + padding);
  }
  std::vector<glm::ivec2> AttemptPacking(const std::vector<glm::ivec2>& rectSizes, int padding, glm::ivec2 size)
  {
    std::vector<glm::ivec2> rectPositions;

    for(size_t rectIndex = 0; rectIndex < rectSizes.size(); rectIndex++)
    {
      bool found = false;
      glm::ivec2 testPos;
      for (testPos.x = 0; testPos.x < size.x && !found; testPos.x++)
      {
        for (testPos.y = 0; testPos.y < size.y && !found; testPos.y++)
        {
          bool ok = (testPos.x >= padding && testPos.x + rectSizes[rectIndex].x <= size.x - padding && testPos.y >= padding && testPos.y + rectSizes[rectIndex].y <= size.y - padding);
          for (size_t testRectIndex = 0; (testRectIndex < rectIndex) && ok; testRectIndex++)
          {
            ok = (ok && TestRects(testPos, rectSizes[rectIndex], rectPositions[testRectIndex], rectSizes[testRectIndex], padding));
          }
          if (ok)
          {
            rectPositions.push_back(testPos);
            found = true;
          }
        }
      }
      if (!found)
        return rectPositions;
    }

    return rectPositions;
  }

  PackingResult PackRects(const std::vector<glm::ivec2>& rectSizes, int padding, int sizeIncrement)
  {
    struct IndexedRect
    {
      glm::ivec2 size;
      size_t index;
    };
    std::vector<IndexedRect> indexedRects;
    for (size_t rectIndex = 0; rectIndex < rectSizes.size(); rectIndex++)
    {
      indexedRects.push_back({rectSizes[rectIndex], rectIndex});
    }
    std::sort(indexedRects.begin(), indexedRects.end(), [](IndexedRect r0, IndexedRect r1) {return r0.size.x * r0.size.y > r1.size.x* r1.size.y; });

    std::vector<glm::ivec2> sortedRects;
    for (auto indexedRect : indexedRects)
    {
      sortedRects.push_back(indexedRect.size);
    }
    glm::ivec2 currSize = glm::ivec2(sizeIncrement, sizeIncrement);

    std::vector<glm::ivec2> sortedPositions;
    while (true)
    {
      sortedPositions = AttemptPacking(sortedRects, padding, currSize);
      if (sortedPositions.size() == rectSizes.size())
        break;
      currSize += glm::ivec2(sizeIncrement, sizeIncrement);
    }

    PackingResult res;
    res.size = currSize;
    res.rectPositions.resize(rectSizes.size());
    for (size_t rectIndex = 0; rectIndex < rectSizes.size(); rectIndex++)
    {

      res.rectPositions[indexedRects[rectIndex].index] = sortedPositions[rectIndex];
    }
    return res;
  }
}