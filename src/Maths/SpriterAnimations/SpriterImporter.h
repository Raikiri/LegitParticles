#pragma once

#include "json/json.h"
#include "../../Utils/GlmInclude.h"

#include <iostream>

namespace almost
{
  namespace spriter
  {
    struct Object
    {
      using FileId = size_t;
      using FolderId = size_t;
      using TimelineId = size_t;
      using TimelineKeyId = size_t;
      using BoneId = size_t;
      using ObjectId = size_t;
      using AnimationId = size_t;

      using MainlineKeyId = size_t;


      struct Animation
      {
        std::string name;
        float length;

        struct TimelineKey
        {
          glm::vec2 relPos = { 0.0f, 0.0f };
          float relAngle = 0;
          float time = 0;
          glm::vec2 relScale = { 0.0f, 0.0f };

          FileId fileId = -1;
          FolderId folderId = -1;

          bool isObject = false;
          bool isBone = false;
        };

        struct MainlineKey
        {
          struct BoneRef
          {
            TimelineId timelineId = -1;
            TimelineKeyId timelineKeyId = -1;
            BoneId parentBoneId = -1;
          };
          std::map<BoneId, BoneRef> boneRefs;

          struct ObjectRef
          {
            BoneId parentBoneId;
            TimelineId timelineId;
            TimelineKeyId timelineKeyId;
            float zIndex;
          };
          std::map<ObjectId, ObjectRef> objectRefs;

          float time;
        };


        struct Timeline
        {
          std::string timelineName;
          std::map<TimelineKeyId, TimelineKey> timelineKeys;
        };
        std::map<TimelineId, Timeline> timelines;

        struct Mainline
        {
          std::map<MainlineKeyId, MainlineKey> mainlineKeys;
        };
        Mainline mainline;
      };
      std::map<AnimationId, Animation> animations;

      struct TextureKey
      {
        FolderId folderId;
        FileId fileId;
        bool operator< (const TextureKey& other) const
        {
          return std::tie(folderId, fileId) < std::tie(other.folderId, other.fileId);
        }
      };
      struct Texture
      {
        std::string filename;
        glm::vec2 pivot;
        glm::vec2 size;
      };
      std::map<TextureKey, Texture> textures;

      static Animation::Timeline ParseTimeline(Json::Value timelineValue);
      static Object::Animation::Mainline ParseMainline(Json::Value mainlineValue);
      static Animation ParseAnimation(Json::Value animationValue);
      static std::map<TextureKey, Texture> ParseTextures(Json::Value folderValue);

      Object() {}
      Object(Json::Value entityValue);
    };

    Object::Animation::Timeline Object::ParseTimeline(Json::Value timelineValue)
    {
      Animation::Timeline timeline;
      timeline.timelineName = timelineValue["name"].asString();

      Json::Value keysArray = timelineValue["key"];
      for (Json::ArrayIndex keyIndex = 0; keyIndex < keysArray.size(); keyIndex++)
      {
        Json::Value keyValue = keysArray[keyIndex];

        Animation::TimelineKey timelineKey;
        Json::Value boneValue = keyValue["bone"];
        Json::Value objectValue = keyValue["object"];

        if (boneValue != Json::Value())
        {
          timelineKey.isBone = true;
        }
        if (objectValue != Json::Value())
        {
          assert(!timelineKey.isBone);
          timelineKey.isObject = true;
          boneValue = objectValue;
        }
        timelineKey.relAngle = boneValue.get("angle", 0.0f).asFloat() / 360.0f * 2.0f * glm::pi<float>();
        timelineKey.relScale = glm::vec2(boneValue.get("scale_x", 1.0f).asFloat(), boneValue.get("scale_y", 1.0f).asFloat());
        timelineKey.relPos = glm::vec2(boneValue.get("x", 0.0f).asFloat(), boneValue.get("y", 0.0f).asInt());
        timelineKey.time = keyValue.get("time", 0.0f).asFloat();
        timelineKey.fileId = FileId(boneValue.get("file", -1).asInt());
        timelineKey.folderId = FolderId(boneValue.get("folder", -1).asInt());

        TimelineKeyId keyId = TimelineKeyId(keyValue["id"].asInt());
        timeline.timelineKeys[keyId] = std::move(timelineKey);
      }
      return timeline;
    }

    Object::Animation::Mainline Object::ParseMainline(Json::Value mainlineValue)
    {
      Animation::Mainline mainline;
      Json::Value keysArray = mainlineValue["key"];
      for (Json::Value keyValue : keysArray)
      {
        Animation::MainlineKey mainlineKey;

        Json::Value boneRefsArray = keyValue["bone_ref"];
        for (Json::Value boneValue : boneRefsArray)
        {
          Animation::MainlineKey::BoneRef boneRef;
          boneRef.parentBoneId = BoneId(boneValue.get("parent", -1).asInt());
          boneRef.timelineId = TimelineId(boneValue["timeline"].asInt());
          boneRef.timelineKeyId = TimelineKeyId(boneValue["key"].asInt());
          BoneId boneId = BoneId(boneValue["id"].asInt());
          mainlineKey.boneRefs[boneId] = boneRef;
        }

        Json::Value objectsArray = keyValue["object_ref"];
        for (Json::Value objectValue : objectsArray)
        {
          Animation::MainlineKey::ObjectRef objectRef;
          objectRef.parentBoneId = BoneId(objectValue.get("parent", -1).asInt());
          objectRef.timelineId = TimelineId(std::stoi(objectValue["timeline"].asString()));
          objectRef.timelineKeyId = TimelineKeyId(objectValue["key"].asInt());
          objectRef.zIndex = objectValue["zIndex"].asFloat();

          ObjectId objectId = ObjectId(objectValue["id"].asInt());
          mainlineKey.objectRefs[objectId] = objectRef;
        }
        mainlineKey.time = keyValue["time"].asFloat();

        MainlineKeyId keyId = MainlineKeyId(keyValue["id"].asInt());
        mainline.mainlineKeys[keyId] = std::move(mainlineKey);
      }
      return mainline;
    }

    Object::Animation Object::ParseAnimation(Json::Value animationValue)
    {
      Animation animation;
      animation.name = animationValue["name"].asString();
      animation.length = animationValue["length"].asFloat();

      std::cout << "Animation name : " << animationValue["name"] << "\n";

      Json::Value timelinesArray = animationValue["timeline"];
      std::cout << "Timelines found : " << timelinesArray.size() << "\n";

      for (Json::ArrayIndex timelineIndex = 0; timelineIndex < timelinesArray.size(); timelineIndex++)
      {
        Json::Value timelineValue = timelinesArray[timelineIndex];
        TimelineId timelineId = TimelineId(timelineValue["id"].asInt());
        animation.timelines[timelineId] = std::move(ParseTimeline(timelineValue));
      }
      animation.mainline = std::move(ParseMainline(animationValue["mainline"]));

      for (const auto &mainlineKey : animation.mainline.mainlineKeys)
      {
        for (const auto &boneRef : mainlineKey.second.boneRefs)
        {
          auto timeline = animation.timelines.find(boneRef.second.timelineId);
          assert(timeline != animation.timelines.end());

          auto timelineKey = timeline->second.timelineKeys.find(boneRef.second.timelineKeyId);
          assert(timelineKey != timeline->second.timelineKeys.end());

          assert(timelineKey->second.isBone);
          //assert(abs(timelineKey->second.time - mainlineKey.second.time) < 1e-2f);
        }

        for (const auto& objectRef : mainlineKey.second.objectRefs)
        {
          auto timeline = animation.timelines.find(objectRef.second.timelineId);
          assert(timeline != animation.timelines.end());
          auto timelineKey = timeline->second.timelineKeys.find(objectRef.second.timelineKeyId);
          assert(timelineKey != timeline->second.timelineKeys.end());

          assert(timelineKey->second.isObject);
          //assert(abs(timelineKey->second.time - mainlineKey.second.time) < 1e-2f);
        }
      }

      return animation;
    }


    std::map<Object::TextureKey, Object::Texture> Object::ParseTextures(Json::Value foldersArray)
    {
      std::map<Object::TextureKey, Object::Texture> textures;
      for (auto folderValue : foldersArray)
      {
        Object::TextureKey textureKey;
        textureKey.folderId = Object::FolderId(folderValue["id"].asInt());
        for (auto fileValue : folderValue["file"])
        {
          textureKey.fileId = Object::FileId(fileValue["id"].asInt());

          Object::Texture tex;
          tex.filename = fileValue["name"].asString();
          tex.pivot = glm::vec2(fileValue.get("pivot_x", 0.5f).asFloat(), fileValue.get("pivot_y", 0.5f).asFloat());
          tex.size = glm::vec2(fileValue.get("width", 1.0f).asFloat(), fileValue.get("height", 1.0f).asFloat());

          textures[textureKey] = tex;
        }
      }
      return textures;
    }


    Object::Object(Json::Value rootValue)
    {
      Json::Value entityValue = rootValue["entity"][0];
      Json::Value animationsArray = entityValue["animation"];

      std::cout << "Animations found : " << animationsArray.size() << "\n";

      for (Json::ArrayIndex animationIndex = 0; animationIndex < animationsArray.size(); animationIndex++)
      {
        Json::Value animationValue = animationsArray[animationIndex];
        AnimationId animationId = AnimationId(animationValue["id"].asInt());
        animations[animationId] = std::move(ParseAnimation(animationValue));
      }

      textures = std::move(ParseTextures(rootValue["folder"]));
    }
  }




  struct FlattenedAnimation
  {
    struct BoneTransform
    {
      glm::vec2 worldPos;
      float worldAngle;
      glm::vec2 worldScale;

      spriter::Object::TextureKey texKey;
    };

    struct Key
    {
      std::vector<BoneTransform> boneTransforms;
    };

    std::vector<Key> keys;
  };

  FlattenedAnimation::BoneTransform LocalTransformToWorld(FlattenedAnimation::BoneTransform parentTransform, FlattenedAnimation::BoneTransform childTransform)
  {
    glm::vec2 xVector = glm::vec2(cos(parentTransform.worldAngle), sin(parentTransform.worldAngle));
    glm::vec2 yVector = glm::vec2(-xVector.y, xVector.x);
    FlattenedAnimation::BoneTransform worldTransform;
    worldTransform.worldAngle = parentTransform.worldAngle + childTransform.worldAngle;
    worldTransform.worldPos = parentTransform.worldPos + xVector * parentTransform.worldScale.x * childTransform.worldPos.x + yVector * parentTransform.worldScale.y * childTransform.worldPos.y;
    worldTransform.worldScale = parentTransform.worldScale * childTransform.worldScale;
    worldTransform.texKey = childTransform.texKey;
    return worldTransform;
  }

  FlattenedAnimation::BoneTransform InterpolateTimeline(const spriter::Object::Animation& anim, spriter::Object::TimelineId timelineId, spriter::Object::TimelineKeyId timelineKeyId, float currTime)
  {

    auto timeline = anim.timelines.find(timelineId);
    assert(timeline != anim.timelines.end());

    auto prevTimelineKey = timeline->second.timelineKeys.find(timelineKeyId);
    assert(prevTimelineKey != timeline->second.timelineKeys.end());

    FlattenedAnimation::BoneTransform prevTransform;
    prevTransform.worldAngle = prevTimelineKey->second.relAngle;
    prevTransform.worldPos = prevTimelineKey->second.relPos;
    prevTransform.worldScale = prevTimelineKey->second.relScale;
    prevTransform.texKey.fileId = prevTimelineKey->second.fileId;
    prevTransform.texKey.folderId = prevTimelineKey->second.folderId;

    float timeWrap = 0.0f;
    auto nextTimelineKey = timeline->second.timelineKeys.find(timelineKeyId + 1);
    if (nextTimelineKey == timeline->second.timelineKeys.end())
    {
      timeWrap = anim.length;
      nextTimelineKey = timeline->second.timelineKeys.find(0);
    }

    if (nextTimelineKey == timeline->second.timelineKeys.end())
      return prevTransform;

    FlattenedAnimation::BoneTransform nextTransform;
    nextTransform.worldAngle = nextTimelineKey->second.relAngle;
    nextTransform.worldPos = nextTimelineKey->second.relPos;
    nextTransform.worldScale = nextTimelineKey->second.relScale;

    float ratio = (currTime - prevTimelineKey->second.time) / (timeWrap + nextTimelineKey->second.time - prevTimelineKey->second.time);
    ratio = glm::clamp(ratio, 0.0f, 1.0f);

    FlattenedAnimation::BoneTransform transform;
    transform.worldAngle = glm::mix(prevTransform.worldAngle, nextTransform.worldAngle, ratio);
    transform.worldPos = glm::mix(prevTransform.worldPos, nextTransform.worldPos, ratio);
    transform.worldScale = glm::mix(prevTransform.worldScale, nextTransform.worldScale, ratio);
    transform.texKey.fileId = prevTransform.texKey.fileId;
    transform.texKey.folderId = prevTransform.texKey.folderId;

    return transform;
  }
  FlattenedAnimation::BoneTransform GetBoneTransform(const spriter::Object::Animation &anim, const spriter::Object::Animation::MainlineKey &mainlineKey, spriter::Object::BoneId boneId, float currTime)
  {
    auto currBone = mainlineKey.boneRefs.find(boneId);
    assert(currBone != mainlineKey.boneRefs.end());

    auto timelineId = currBone->second.timelineId;
    auto timelineKeyId = currBone->second.timelineKeyId;

    auto currTransform = InterpolateTimeline(anim, timelineId, timelineKeyId, currTime);
       //= GetBoneTransform(anim, mainlineKey, boneId, currTime);
    if (currBone->second.parentBoneId != spriter::Object::BoneId(-1))
    {
      FlattenedAnimation::BoneTransform parentTransform = GetBoneTransform(anim, mainlineKey, currBone->second.parentBoneId, currTime);
      return LocalTransformToWorld(parentTransform, currTransform);
    }
    return currTransform;
  }

  FlattenedAnimation ConvertToFlattenedAnimation(const spriter::Object::Animation &anim, float timeStep)
  {
    auto currMainlineKey = anim.mainline.mainlineKeys.begin();
    FlattenedAnimation flatAnimation;

    /*for (auto objRef : currMainlineKey->second.objectRefs)
    {
      auto timelineId = objRef.second.timelineId;
      auto timelineKeyId = objRef.second.timelineKeyId;

      auto timeline = anim.timelines.find(timelineId);
      assert(timeline != anim.timelines.end());

      auto prevTimelineKey = timeline->second.timelineKeys.find(timelineKeyId);
      assert(prevTimelineKey != timeline->second.timelineKeys.end());

      flatAnimation.bones.
    }*/


    for (float currTime = 0.0f; currTime < anim.length; currTime += timeStep)
    {
      FlattenedAnimation::Key flattenedKey;
      auto nextMainlineKey = std::next(currMainlineKey);
      while (nextMainlineKey != anim.mainline.mainlineKeys.end() && nextMainlineKey->second.time < currTime)
      {
        currMainlineKey = nextMainlineKey;
        nextMainlineKey = std::next(currMainlineKey);
      }

      for (auto objRef : currMainlineKey->second.objectRefs)
      {
        auto objTransform = InterpolateTimeline(anim, objRef.second.timelineId, objRef.second.timelineKeyId, currTime);
        auto parentTransform = GetBoneTransform(anim, currMainlineKey->second, objRef.second.parentBoneId, currTime);

        objTransform.worldAngle = 0.0f;
        objTransform.worldPos = glm::vec2(0.0f, 0.0f);
        flattenedKey.boneTransforms.push_back(LocalTransformToWorld(parentTransform, objTransform));
      }
      flatAnimation.keys.push_back(flattenedKey);
    }
    return flatAnimation;
  }

}