@echo off
setlocal EnableExtensions EnableDelayedExpansion
set CompilerExe="d:\Portable soft\VulkanSDK\1.1.121.2\Bin\glslangValidator.exe"
set OptimizerExe="d:\Portable soft\VulkanSDK\1.1.121.2\Bin\spirv-opt.exe"
set OptimizerConfig="OptimizerConfig.cfg"

set inname=d:\CPP Projects\LegitEngine\bin\data\Shaders\glsl\PointRenderer\listBlockSizedGathering.frag
rem set inname=d:\CPP Projects\LegitEngine\bin\data\Shaders\glsl\PointRenderer\listBlockSizedGathering - arg.frag
rem set inname=d:\CPP Projects\LegitEngine\bin\data\Shaders\glsl\PointRenderer\arg test.frag
set outname=!inname:\glsl\=\spirv\!
@echo Buiding !inname! to !outname!
%CompilerExe% -V "!inname!" -l --target-env vulkan1.1 -o "!outname!".spv
rem %OptimizerExe% "!outname!"_u.spv -Oconfig="OptimizerConfig.cfg" -o "!outname!".spv
