#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(binding = 0, set = 0) uniform GlobalDataBuffer
{
	vec2 camPos;
  vec2 camFov;
  vec2 viewportSize;
  float camAngle;
  float debugVal;
}globalDataBuf;

/*layout(binding = 0, set = 1) uniform DrawCallData
{
	mat4 modelMatrix; //object->world
	vec4 albedoColor;
	vec4 emissiveColor;
};*/
uniform layout(binding = 1, set = 0) sampler2D animationTex;

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
  int rowsPerBone = 5;
  int keyIndex = int(globalDataBuf.debugVal * 500.0f);

  outColor.rgba = vec4(0.0f);
  for(int boneIndex = 0; boneIndex < 15; boneIndex++)
  {
    //int boneIndex = 7;
    vec4 posScaleSample = texelFetch(animationTex, ivec2(keyIndex, boneIndex * rowsPerBone + 0), 0);
    vec4 angleSample = texelFetch(animationTex, ivec2(keyIndex, boneIndex * rowsPerBone + 1), 0);
    vec4 uvSample = texelFetch(animationTex, ivec2(keyIndex, boneIndex * rowsPerBone + 2), 0);
    
    vec2 uvMin = uvSample.xy;
    vec2 uvMax = uvSample.zw;
    
    vec2 pos = (posScaleSample.xy - 0.5f) * 400.0f;
    vec2 scale = posScaleSample.zw * 400.0f * 2.0f;
    
    float ang = angleSample.r * 6.28f;
    
    vec2 xVec = vec2(cos(ang), sin(ang));
    vec2 yVec = vec2(-xVec.y, xVec.x);
    
    vec2 pixelPos = gl_FragCoord.xy - globalDataBuf.viewportSize.xy / 2;
    
    vec2 localUv = vec2(dot(xVec, pixelPos - pos) / scale.x, dot(yVec, pixelPos - pos) / scale.y) + vec2(0.25f);

    if(localUv.x > 0.0f && localUv.y > 0.0f && localUv.x < 1.0f && localUv.y < 1.0f)
    {
      vec2 uv = uvMin + (uvMax - uvMin) * localUv;
      
      vec4 layerColor = texture(animationTex, uv);
      outColor.rgba = mix(outColor.rgba, layerColor, layerColor.a);
      //outColor.rgba = vec4(localUv, 0.0f, 1.0f);
    }
  }
  /*vec2 screenPos = gl_FragCoord.xy / globalDataBuf.viewportSize.xy;
  vec2 uv = uvMin + (uvMax - uvMin) * screenPos;
  outColor.rgba = texture(animationTex, uv);*/
  return;
}