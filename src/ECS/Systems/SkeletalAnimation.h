#pragma once
#include "../../LegitVulkan/ProfilerTask.h"
#include "../../LegitVulkan/CpuProfiler.h"

#include <entity\registry.hpp>
#include "../../Utils/GroupArg.h"
#include "../Context/AnimationData.h"
#include "json/json.h"
#include "../../Maths/RectPacker/RectPacker.h"
#include "lodepng.h"

namespace almost
{
  //const glm::vec4 favBlue = glm::vec4(glm::pow(glm::vec3(0.0, 0.5, 0.75), glm::vec3(2.2f)), 1.0f);
  //const glm::vec4 favOrange = glm::vec4(glm::pow(glm::vec3(0.8f, 0.4f, 0.1f), glm::vec3(2.2f)), 1.0f);

  struct PackedTexture
  {
    glm::ivec2 pos;
    glm::ivec2 size;
    std::string filename;
  };

  unsigned char EncodeUNorm(float val)
  {
    return (unsigned char)(glm::clamp<int>(int(255 * val), 0, 255));
  }

  void PackTextures(std::vector<PackedTexture> packedTextures, glm::ivec2 dstSize, FlattenedAnimation animation, std::map<spriter::Object::TextureKey, size_t> textureIndices, std::string srcPrefix, std::string dstFilename)
  {
    int rowsPerBone = 5;
    int metadataSize = int(packedTextures.size() * rowsPerBone);

    glm::ivec2 finalSize = dstSize + glm::ivec2(0, metadataSize);

    std::vector<unsigned char> dstData;
    dstData.resize(sizeof(glm::uint32_t) * finalSize.x * finalSize.y);
    for (int y = 0; y < finalSize.y; y++)
    {
      for (int x = 0; x < finalSize.x; x++)
      {
        for (int channel = 0; channel < 4; channel++)
          dstData[(x + y * finalSize.x) * 4 + channel] = 0;
      }
    }

    size_t texIndex = 0;
    for (size_t keyIndex = 0; keyIndex < std::min(animation.keys.size(), size_t(dstSize.x)); keyIndex++)
    {
      auto key = animation.keys[keyIndex];
      for (size_t boneIndex = 0; boneIndex < key.boneTransforms.size(); boneIndex++)
      {
        auto bone = key.boneTransforms[boneIndex];
        auto texKey = bone.texKey;
        auto texIndexIt = textureIndices.find(texKey);
        assert(texIndexIt != textureIndices.end());

        auto packedTex = packedTextures[texIndexIt->second];


        dstData[(keyIndex + (boneIndex * rowsPerBone + 0) * finalSize.x) * 4 + 0] = EncodeUNorm(bone.worldPos.x / 400.0f + 0.5f);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 0) * finalSize.x) * 4 + 1] = EncodeUNorm(bone.worldPos.y / 400.0f + 0.5f);
        glm::vec2 pixelScale = bone.worldScale * glm::vec2(packedTex.size);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 0) * finalSize.x) * 4 + 2] = EncodeUNorm(pixelScale.x / 400.0f);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 0) * finalSize.x) * 4 + 3] = EncodeUNorm(pixelScale.y / 400.0f);

        dstData[(keyIndex + (boneIndex * rowsPerBone + 1) * finalSize.x) * 4 + 0] = EncodeUNorm(glm::fract(bone.worldAngle / 6.28f));
        dstData[(keyIndex + (boneIndex * rowsPerBone + 1) * finalSize.x) * 4 + 1] = EncodeUNorm(0.0f);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 1) * finalSize.x) * 4 + 2] = EncodeUNorm(0.0f);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 1) * finalSize.x) * 4 + 3] = EncodeUNorm(0.0f);


        glm::vec2 texPosMin = glm::vec2(glm::ivec2(0, metadataSize) + packedTex.pos) / glm::vec2(finalSize);
        glm::vec2 texPosMax = glm::vec2(glm::ivec2(0, metadataSize) + packedTex.pos + packedTex.size) / glm::vec2(finalSize);

        dstData[(keyIndex + (boneIndex * rowsPerBone + 2) * finalSize.x) * 4 + 0] = EncodeUNorm(texPosMin.x);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 2) * finalSize.x) * 4 + 1] = EncodeUNorm(texPosMin.y);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 2) * finalSize.x) * 4 + 2] = EncodeUNorm(texPosMax.x);
        dstData[(keyIndex + (boneIndex * rowsPerBone + 2) * finalSize.x) * 4 + 3] = EncodeUNorm(texPosMax.y);

      }
    }
    for (auto packedTex : packedTextures)
    {
      std::vector<unsigned char> srcData;
      glm::uvec2 srcSize;
      auto res = lodepng::decode(srcData, srcSize.x, srcSize.y, srcPrefix + packedTex.filename);
      assert(!res);
      assert(packedTex.size.x == srcSize.x && packedTex.size.y == srcSize.y);


      glm::ivec2 srcPixel;
      for (srcPixel.y = 0; srcPixel.y < packedTex.size.y; srcPixel.y++)
      {
        for (srcPixel.x = 0; srcPixel.x < packedTex.size.x; srcPixel.x++)
        {
          glm::ivec2 dstPixel = srcPixel + packedTex.pos + glm::ivec2(0, metadataSize);

          for (int channel = 0; channel < 4; channel++)
            dstData[(dstPixel.x + dstPixel.y * finalSize.x) * 4 + channel] = srcData[(srcPixel.x + srcPixel.y * srcSize.x) * 4 + channel];
        }
      }
      texIndex++;
    }
    lodepng::encode(dstFilename, (unsigned char *)dstData.data(), finalSize.x, finalSize.y);
  }

  AnimationData InitAnimationData()
  {
    AnimationData animationData;

    std::string folder = "../data/Animations/WalkTest/";
    std::string animName = "test.scon";
    std::ifstream fileStream(folder + animName);
    if (!fileStream.is_open())
    {
      std::cout << "Can't open animation file " << animName << "\n";
      return animationData;
    }

    Json::Value animRoot;
    Json::Reader reader;

    std::cout << "Parsing animation file " << animName << "\n";
    bool result = reader.parse(fileStream, animRoot);

    if (result)
    {
      std::cout << "File " << animName << ", parsing successful\n";
    }
    else
    {
      std::cout << "Error: File " << animName << ", parsing failed with errors: " << reader.getFormattedErrorMessages() << "\n";
      return animationData;
    }

    animationData.spriterObject = spriter::Object(animRoot);

    for (auto anim : animationData.spriterObject.animations)
    {
      if (anim.second.name == "walk")
      {
        animationData.flattenedAnimation = ConvertToFlattenedAnimation(anim.second, anim.second.length / 128.0f);
      }
    }


    size_t texIndex = 0;
    std::vector<PackedTexture> packedTextures;
    std::map<spriter::Object::TextureKey, size_t> textureIndices;
    std::vector<glm::ivec2> rectSizes;
    for (auto key : animationData.flattenedAnimation.keys)
    {
      for (auto bone : key.boneTransforms)
      {
        auto texIt = textureIndices.find(bone.texKey);
        if (texIt == textureIndices.end())
        {
          auto baseTexIt = animationData.spriterObject.textures.find(bone.texKey);
          assert(baseTexIt != animationData.spriterObject.textures.end());

          PackedTexture packedTex;
          packedTex.filename = baseTexIt->second.filename;
          packedTex.pos = glm::ivec2(-1, -1);
          packedTex.size = baseTexIt->second.size;
          packedTextures.push_back(packedTex);

          rectSizes.push_back(glm::ivec2(baseTexIt->second.size));
          textureIndices[bone.texKey] = texIndex++;
        }
      }
    }

    auto packingRes = PackRects(rectSizes, 2, 32);
    for (size_t texIndex = 0; texIndex < textureIndices.size(); texIndex++)
    {
      packedTextures[texIndex].pos = packingRes.rectPositions[texIndex];
    }

    PackTextures(packedTextures, packingRes.size, animationData.flattenedAnimation, textureIndices, "../data/Animations/WalkTest/", "../data/Out/atlas.png");
    return animationData;
  }
  struct WindowData;
  struct InputData;
  struct CameraData;
  struct MeshRendererData;

  
  void SubmitAnimation(
    almost::AnimationData& animationData,
    almost::MeshRendererData& meshRendererData)
  {
    animationData.currKeyIndex = ((animationData.currKeyIndex + 1) % animationData.flattenedAnimation.keys.size());
    auto currKey = animationData.flattenedAnimation.keys[animationData.currKeyIndex];
    for (auto boneTransform : currKey.boneTransforms)
    {
      glm::vec2 xVector = glm::vec2(cos(boneTransform.worldAngle), sin(boneTransform.worldAngle));
      glm::vec2 yVector = glm::vec2(-xVector.y, xVector.x);
      std::array<glm::vec2, 3> points;
      glm::vec2 localSize = { 200.0f, 5.0f };
      points[0] = boneTransform.worldPos - yVector * boneTransform.worldScale.y * localSize.y;
      points[1] = boneTransform.worldPos + yVector * boneTransform.worldScale.y * localSize.y;
      points[2] = boneTransform.worldPos + xVector * boneTransform.worldScale.x * localSize.x;
      almost::SubmitTriangle(meshRendererData, points, favBlue);
    }
  }
}