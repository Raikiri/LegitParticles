#include "../Context/CameraData.h"
#include "../Context/RendererData.h"
#include "../Context/WindowData.h"
#include "../Context/InputData.h"

namespace almost
{

  CameraData InitCameraData()
  {
    CameraData cameraData;
    cameraData.pos = glm::vec2(0.0f, 0.0f);
    cameraData.ang = 0.0f;
    cameraData.fov = glm::vec2(500.0f, 1.0f);
    //cameraData.fov = glm::
    return cameraData;
  }

  void UpdateCamera(const WindowData &windowData, const RendererData &rendererData, InputData &inputData, CameraData& cameraData)
  {
    glm::vec2 cameraRight = glm::vec2(cos(cameraData.ang), sin(cameraData.ang));
    glm::vec2 cameraUp = glm::vec2(-cameraRight.y, cameraRight.x);

    glm::vec2 dir = glm::vec2(0.0f);
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_E))
      dir += cameraUp;
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_S))
      dir += -cameraRight;
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_D))
      dir += -cameraUp;
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_F))
      dir += cameraRight;

    float angSpeed = 1.0f;
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_R))
      cameraData.ang -= angSpeed * inputData.deltaTime;
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_W))
      cameraData.ang += angSpeed * inputData.deltaTime;


    float scaleSpeed = 1.0f;
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_T))
      cameraData.fov.x *= exp(-inputData.deltaTime * scaleSpeed);
    if (glfwGetKey(windowData.window->glfw_window, GLFW_KEY_G))
      cameraData.fov.x *= exp(inputData.deltaTime * scaleSpeed);

    float moveSpeed = cameraData.fov.x * 1.0f;
    cameraData.pos += dir * moveSpeed * inputData.deltaTime;

    if (rendererData.inFlightQueue)
    {
      auto extent = rendererData.inFlightQueue->GetImageSize();
      cameraData.fov = GetFov2(cameraData.fov.x, glm::vec2(extent.width, extent.height));

      glm::vec2 screenPos = glm::vec2(inputData.mousePos) / glm::vec2(extent.width, extent.height);
      screenPos.y = 1.0f - screenPos.y;
      inputData.prevWorldMousePos = inputData.worldMousePos;
      inputData.worldMousePos = Unproject(cameraData.pos, cameraData.ang, cameraData.fov, screenPos);
    }

    //cameraData.fov.y = rendererData.
  }

  void ProcessCameraInput()
  {
    /*if (!ImGui::IsWindowFocused(ImGuiFocusedFlags_AnyWindow))
    {
      glm::vec3 dir = glm::vec3(0.0f, 0.0f, 0.0f);

      if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2))
      {
        float mouseSpeed = 0.01f;
        camera.horAngle += float((mousePos - prevMousePos).x * mouseSpeed);
        camera.vertAngle += float((mousePos - prevMousePos).y * mouseSpeed);
      }
      glm::mat4 cameraTransform = camera.GetTransformMatrix();
      glm::vec3 cameraForward = glm::vec3(cameraTransform * glm::vec4(0.0f, 0.0f, 1.0f, 0.0f));
      glm::vec3 cameraRight = glm::vec3(cameraTransform * glm::vec4(1.0f, 0.0f, 0.0f, 0.0f));
      glm::vec3 cameraUp = glm::vec3(cameraTransform * glm::vec4(0.0f, 1.0f, 0.0f, 0.0f));


      if (glfwGetKey(window, GLFW_KEY_E))
        dir += glm::vec3(0.0f, -0.0f, 1.0f);
      if (glfwGetKey(window, GLFW_KEY_S))
        dir += glm::vec3(-1.0f, 0.0f, 0.0f);
      if (glfwGetKey(window, GLFW_KEY_D))
        dir += glm::vec3(0.0f, 0.0f, -1.0f);
      if (glfwGetKey(window, GLFW_KEY_F))
        dir += glm::vec3(1.0f, 0.0f, 0.0f);
      if (glfwGetKey(window, GLFW_KEY_SPACE))
        dir += glm::vec3(0.0f, 1.0f, 0.0f);
      if (glfwGetKey(window, GLFW_KEY_C))
        dir += glm::vec3(0.0f, -1.0f, 0.0f);

      float cameraSpeed = 3.0f;
      camera.pos += cameraForward * dir.z * cameraSpeed * deltaTime;
      camera.pos += cameraRight * dir.x * cameraSpeed * deltaTime;
      camera.pos += cameraUp * dir.y * cameraSpeed * deltaTime;
    }*/
  }
}
