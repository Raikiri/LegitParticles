#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(location = 0) in vec3 attribPosition;
layout(location = 1) in vec2 attribUv;
layout(location = 2) in vec4 attribColor;

layout(binding = 0, set = 0) uniform GlobalDataBuffer
{
	mat4 viewProjMatrix; //world->screen
} globalDataBuf;

/*layout(binding = 0, set = 1) uniform DrawCallData
{
	mat4 modelMatrix; //object->world
	vec4 albedoColor;
	vec4 emissiveColor;
};*/

out gl_PerVertex 
{
	vec4 gl_Position;
};

layout(location = 0) out vec2 vertUv;
layout(location = 1) out vec4 vertColor;

void main()
{
	gl_Position = globalDataBuf.viewProjMatrix * vec4(attribPosition, 1.0f);
	vertUv = attribUv;
	vertColor = attribColor;
}