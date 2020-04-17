#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(binding = 0, set = 0) uniform GlobalDataBuffer
{
	vec2 camPos;
  vec2 camFov;
  vec2 viewportSize;
  float camAngle;
}globalDataBuf;

/*layout(binding = 0, set = 1) uniform DrawCallData
{
	mat4 modelMatrix; //object->world
	vec4 albedoColor;
	vec4 emissiveColor;
};*/

layout(location = 0) in vec2 fragScreenPos;

layout(location = 0) out vec4 outColor;
float saturate(float val)
{
  return clamp(val, 0.0f, 1.0f);
}
float GetGridDist(vec2 pos, float step)
{
  vec2 phases = fract(pos / step + vec2(0.5f));
  vec2 gridDists = 1.0f - min(phases * 2.0f, vec2(2.0f) - phases * 2.0f);
  return min(gridDists.x, gridDists.y) * step * 0.5f;
}

float GetGrid(vec2 pos, float step, float width, float aaWidth)
{
  float gridDist = GetGridDist(pos, step);
  return saturate(1.0f - ((gridDist - width * 0.5) / aaWidth + 0.5));
}


void main() 
{
  vec2 camRightVec = vec2(cos(globalDataBuf.camAngle), sin(globalDataBuf.camAngle));
  vec2 camUpVec = vec2(-camRightVec.y, camRightVec.x);
  vec2 projectedCoord = (fragScreenPos * 2.0f - vec2(1.0f));
  vec2 worldPos = globalDataBuf.camPos + camRightVec * globalDataBuf.camFov.x * projectedCoord.x + camUpVec * globalDataBuf.camFov.y * projectedCoord.y;
  
  
  float aaWidth = 3.0f * globalDataBuf.camFov.x / globalDataBuf.viewportSize.x;
  float totalGridIntensity = 0.0f;
  
  float gridMult = 10.0f;
  const int gridsCount = 2;
  float targetMaxStep = globalDataBuf.camFov.x;
  float targetMinStep = targetMaxStep / pow(gridMult, gridsCount);
  int maxPow = int(floor(log(targetMaxStep) / log(gridMult)));
  int minPow = maxPow - 3;
  for(int gridIndex = 0; gridIndex < gridsCount; gridIndex++)
  {
    float gridStep = pow(gridMult, maxPow - gridIndex);
    float width = gridStep * 0.02f;

    float gridOpacity = saturate((1.0f - gridStep / targetMaxStep) * 5.0f);
    gridOpacity *= saturate((gridStep / targetMinStep - 1.0f) * 5.0f);
    float gridIntensity = GetGrid(worldPos, gridStep, width, aaWidth) * gridOpacity;
    totalGridIntensity = max(gridIntensity, totalGridIntensity);
  }
  vec4 gridBackground = vec4(vec3(0.03f), 1.0f);
  vec4 gridForeground = vec4(vec3(0.07f), 1.0f);
  outColor = mix(gridBackground, gridForeground, totalGridIntensity);
}