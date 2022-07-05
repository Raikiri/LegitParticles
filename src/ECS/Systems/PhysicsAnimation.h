#pragma once

#include <entity\registry.hpp>
#include "../../Utils/GroupArg.h"
#include "../Context/PhysicsAnimationData.h"
#include "json/json.h"

namespace almost
{

  glm::vec2 GetNextPhysicsAnimationPosNumerical(glm::vec2 prevAnimPos, glm::vec2 animPos, glm::vec2 prevControlPos, glm::vec2 controlPos, float freq, float damp, float response, float dt)
  {
    glm::vec2 animVelocity = (animPos - prevAnimPos) / dt;
    glm::vec2 controlVelocity = (controlPos - prevControlPos) / dt;
    float pi = 3.141592f;
    float omega = 2.0f * pi * freq;
    glm::vec2 accel = (controlPos - animPos) * omega * omega;
    accel += response * damp * omega * controlVelocity - 2.0f * omega * damp * animVelocity;

    return animPos + (animPos - prevAnimPos) + accel * dt * dt;
  }

  template<typename Vec>
  using ComplexVec = glm::vec<2, Vec>;

  template<typename Vec>
  ComplexVec<Vec> ComplexMul(ComplexVec<Vec> a, ComplexVec<Vec> b)
  {
    return ComplexVec<Vec>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
  }

  template<typename Vec>
  ComplexVec<Vec> ComplexExp(Vec ang)
  {
    return ComplexVec<Vec>(glm::cos(ang), glm::sin(ang));
  }

  template<typename Vec>
  ComplexVec<Vec> ComplexConj(ComplexVec<Vec> a)
  {
    return ComplexVec<Vec>(a.x, -a.y);
  }

  /*SmoothState AdvanceSmoothState_(SmoothState animState, SmoothState controlState, float freq, float damp, float response, float dt)
  {
    using Complex = ComplexVec<glm::vec2>;
    glm::vec2 animVelocity = animState.velocity;
    glm::vec2 controlVelocity = controlState.velocity;
    float pi = 3.141592f;
    float omega = 2.0f * pi * freq;

    glm::vec2 animPosition = animState.position - controlState.position;

    Complex A0 = Complex(animPosition / 2.0f, -animVelocity / omega / 2.0f);
    Complex A1 = ComplexConj(A0);

    //amplitude = -animVelocity / omega / glm::sin(phi);

    Complex phaseAdvance0 = ComplexExp(glm::vec2( omega * dt));
    Complex phaseAdvance1 = ComplexConj(phaseAdvance0);

    Complex velocityMult0 = Complex(glm::vec2(0.0f),  glm::vec2(omega));
    Complex velocityMult1 = ComplexConj(velocityMult0);

    Complex advancedPos = ComplexMul(A0, phaseAdvance0) + ComplexMul(A1, phaseAdvance1);
    Complex advancedVelocity = ComplexMul(velocityMult0, ComplexMul(A0, phaseAdvance0)) + ComplexMul(velocityMult1, ComplexMul(A1, phaseAdvance1));

    SmoothState advancedState;
    advancedState.position = controlState.position + advancedPos.x;
    advancedState.velocity = advancedVelocity.x;
    return advancedState;
  }*/

  SmoothState AdvanceSmoothState(SmoothState animState, SmoothState controlState, float freq, float damp, float response, float dt)
  {
    damp = glm::min(0.999f, damp);
    using Complex = ComplexVec<glm::vec2>;
    glm::vec2 animVelocity = animState.velocity;
    glm::vec2 controlVelocity = controlState.velocity;
    float pi = 3.141592f;
    float omega = 2.0f * pi * freq;

    glm::vec2 deltaControlPos = response * damp * controlVelocity / omega;
    glm::vec2 animPosition = animState.position - (controlState.position + deltaControlPos);

    float s = glm::sqrt(glm::max(0.0f, 1.0f - damp * damp));
    Complex A0 = Complex(animPosition / 2.0f, -(animVelocity + animPosition * damp * omega) / (omega * 2.0f * s));
    Complex A1 = ComplexConj(A0);

    //amplitude = -animVelocity / omega / glm::sin(phi);

    Complex phaseAdvance0 = ComplexExp(glm::vec2(omega * s * dt));
    Complex phaseAdvance1 = ComplexConj(phaseAdvance0);

    Complex velocityMult0 = Complex(glm::vec2(0.0f), glm::vec2(omega * s));
    Complex velocityMult1 = ComplexConj(velocityMult0);

    Complex dampExp = glm::vec2(glm::exp(-damp * omega * dt));
    Complex advancedPos = (ComplexMul(A0, phaseAdvance0) + ComplexMul(A1, phaseAdvance1)) * dampExp;
    Complex advancedVelocity = -advancedPos * glm::vec2(damp * omega) + 
      (ComplexMul(velocityMult0, ComplexMul(A0, phaseAdvance0)) + ComplexMul(velocityMult1, ComplexMul(A1, phaseAdvance1))) * dampExp;

    SmoothState advancedState;
    advancedState.position = controlState.position + deltaControlPos + advancedPos.x;
    advancedState.velocity = advancedVelocity.x;
    return advancedState;
  }

  //const glm::vec4 favBlue = glm::vec4(glm::pow(glm::vec3(0.0, 0.5, 0.75), glm::vec3(2.2f)), 1.0f);
  //const glm::vec4 favOrange = glm::vec4(glm::pow(glm::vec3(0.8f, 0.4f, 0.1f), glm::vec3(2.2f)), 1.0f);

  PhysicsAnimationData InitPhysicsAnimationData()
  {
    PhysicsAnimationData physicsAnimationData;
    physicsAnimationData.objPos0 = glm::vec2(0.0f, 0.0f);
    physicsAnimationData.prevObjPos0 = glm::vec2(0.0f, 0.0f);

    physicsAnimationData.objState1.position = glm::vec2(0.0f, 0.0f);
    physicsAnimationData.objState1.velocity = glm::vec2(0.0f, 0.0f);

    physicsAnimationData.prevMouseWorldPos = glm::vec2(0.0f, 0.0f);

    physicsAnimationData.lastUpdateTime = std::chrono::system_clock::now();
    return physicsAnimationData;
  }
  struct WindowData;
  struct InputData;
  struct CameraData;
  struct MeshRendererData;

  void ProcessPhysicsAnimation(
    almost::PhysicsAnimationData& physicsAnimationData,
    almost::InputData& inputData,
    almost::WindowData& windowData)
  {
    static float freq = 1.0f;
    ImGui::SliderFloat("Freq", &freq, 0.00001f, 10.0f);
    static float damp = 1.0f;
    ImGui::SliderFloat("Damp", &damp, 0.0f, 2.0f);
    static float response = 1.0f;
    ImGui::SliderFloat("Response", &response, 0.0f, 2.0f);
    glm::vec2 nextObjPos;
    static bool analytical = false;
    ImGui::Checkbox("Analytical", &analytical);
    //float dt = 1.0f / 240.0f;
    /*if (analytical)
      nextObjPos = GetNextPhysicsAnimationPosAnalytical(physicsAnimationData.prevObjPos, physicsAnimationData.objPos, inputData.prevWorldMousePos, inputData.worldMousePos, freq, damp, response, dt);
    else*/

    /*SmoothState testObjState;
    testObjState.position = glm::vec2(100.0f, 0.0f);
    testObjState.velocity = glm::vec2(0.0f, 0.0f);

    SmoothState testControlState;
    testControlState.position = glm::vec2(0.0f, 0.0f);
    testControlState.velocity = glm::vec2(0.0f, 0.0f);

    for (int i = 0; i < 100; i++)
    {
      testObjState = AdvanceSmoothState(testObjState, testControlState, 1.0f, 1.0f, 0.0f, 0.25f);
      int p = 1;
    }*/



    static int updateTimeMs = 10;
    ImGui::SliderInt("Dt (ms)", &updateTimeMs, 1, 300);

    auto updatePeriod = std::chrono::milliseconds(int(updateTimeMs));
    auto currTime = std::chrono::system_clock::now();
    if (std::chrono::duration_cast<std::chrono::milliseconds>(currTime - physicsAnimationData.lastUpdateTime) > updatePeriod)
    {
      float dt = updateTimeMs / 1000.0f;
      physicsAnimationData.lastUpdateTime += updatePeriod;
      {
        nextObjPos = GetNextPhysicsAnimationPosNumerical(physicsAnimationData.prevObjPos0, physicsAnimationData.objPos0, physicsAnimationData.prevMouseWorldPos, inputData.worldMousePos, freq, damp, response, float(updateTimeMs) / 1000.0f);

        physicsAnimationData.prevObjPos0 = physicsAnimationData.objPos0;
        physicsAnimationData.objPos0 = nextObjPos;
      }

      {
        SmoothState controlState;
        controlState.position = inputData.worldMousePos;
        controlState.velocity = (inputData.worldMousePos - physicsAnimationData.prevMouseWorldPos) / dt;
        physicsAnimationData.objState1 = AdvanceSmoothState(physicsAnimationData.objState1, controlState, freq, damp, response, float(updateTimeMs) / 1000.0f);
      }
      physicsAnimationData.prevMouseWorldPos = inputData.worldMousePos;
    }

    if (glfwGetKey(windowData.window, GLFW_KEY_SPACE))
    {
      physicsAnimationData.objPos0 = inputData.worldMousePos;
      physicsAnimationData.prevObjPos0 = inputData.worldMousePos;
      physicsAnimationData.objState1.position = inputData.worldMousePos;
      physicsAnimationData.objState1.velocity = glm::vec2(0.0f, 0.0f);
    }
  }
  
  void SubmitPhysicsAnimation(
    almost::PhysicsAnimationData& physicsAnimationData,
    almost::InputData& inputData,
    almost::MeshRendererData& meshRendererData)
  {
    almost::SubmitCircle(meshRendererData, physicsAnimationData.objState1.position, 15.0f, 10, favBlue);
    almost::SubmitCircle(meshRendererData, physicsAnimationData.objPos0, 10.0f, 10, favOrange);
  }
}