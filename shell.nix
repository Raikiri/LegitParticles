with import <nixpkgs> {};

mkShell {
  name = "dotnet";
  packages = [
    wayland
    glfw
    cmake
    clang
    lldb
    gnumake
    powershell
    vulkan-headers
    vulkan-loader
    vulkan-validation-layers
    vulkan-tools        # vulkaninfo
    glslang
    shaderc             # GLSL to SPIRV compiler - glslc
    renderdoc           # Graphics debugger
    tracy               # Graphics profiler
    vulkan-tools-lunarg # vkconfig
  ];

  buildInputs = with pkgs; [
  ];


  
  LD_LIBRARY_PATH="${glfw}/lib:${freetype}/lib:${vulkan-loader}/lib:${vulkan-headers}/lib:${vulkan-validation-layers}";
  #VULKAN_SDK_PATH = "${vulkan-headers}";
  
  #this fixes glfw wayland initialization error by using x11 instead
  shellHook = ''
          export LD_LIBRARY_PATH=${pkgs.wayland}/lib:$LD_LIBRARY_PATH
          '';
}


