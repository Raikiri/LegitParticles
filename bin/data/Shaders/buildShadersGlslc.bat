@echo off
setlocal EnableExtensions EnableDelayedExpansion
set CompilerExe="d:\Portable soft\VulkanSDK\1.1.121.2\Bin\glslc.exe"

for /r glsl/ %%I in (*.vert) do (
  set outname=%%I.spv 
  set outname=!outname:\glsl\=\spirv\!
  @echo Building %%I
  @echo To !outname!
  %CompilerExe% "%%I" -o "!outname!"
)
for /r glsl/ %%I in (*.frag) do (
  set outname=%%I.spv
  set outname=!outname:\glsl\=\spirv\!
  @echo Building %%I
  @echo To !outname!
  %CompilerExe% "%%I" -o "!outname!"
)

for /r glsl/ %%I in (*.comp) do (
  set outname=%%I.spv
  set outname=!outname:\glsl\=\spirv\!
  @echo Building %%I
  @echo To !outname!
  %CompilerExe% "%%I" -o "!outname!"
)
rem pause

rem 
rem %CompilerExe% -V glsl/fragmentShader.frag -o spirv/fragmentShader.spv
rem %CompilerExe% -V glsl/frameDescriptorLayout.comp -o spirv/frameDescriptorLayout.spv
rem %CompilerExe% -V glsl/passDescriptorLayout.comp -o spirv/passDescriptorLayout.spv
