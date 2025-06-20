cmake_minimum_required(VERSION 3.29.0)
project(LegitEngineCMake VERSION 0.1.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(Vulkan)
message(STATUS "Found Vulkan in ${Vulkan_INCLUDE_DIRS}")
message(STATUS "Found Vulkan libs in ${Vulkan_LIBRARIES}")

set(TINYOBJ_INCLUDE_DIR ./dependencies/tiny_obj_loader/include)
set(SPIRVCROSS_INCLUDE_DIR ./dependencies/spirv-cross)
set(GLI_INCLUDE_DIR ./dependencies/gli)
set(ENTT_INCLUDE_DIR ./dependencies/entt)
set(GLM_INCLUDE_DIR ./dependencies/glm/include)
set(JSONCPP_INCLUDE_DIR ./dependencies/jsoncpp)
set(IMGUI_INCLUDE_DIR ./dependencies/imgui)

#set(GLFW3_LIBRARIES "${CMAKE_SOURCE_DIR}/dependencies/glfw/lib-vc2015/glfw3.lib")
#find_package(GLM)
#find_package(GLFW3)
set(GLFW_USE_WAYLAND ON)
find_package(GLFW3)

include_directories ("${Vulkan_INCLUDE_DIRS}")
include_directories ("${ENTT_INCLUDE_DIR}")
include_directories ("${GLM_INCLUDE_DIR}")
include_directories ("${GLFW3_INCLUDE_DIR}")
include_directories ("${SPIRVCROSS_INCLUDE_DIR}")
include_directories ("${GLI_INCLUDE_DIR}")
include_directories ("${TINYOBJ_INCLUDE_DIR}")
include_directories ("${JSONCPP_INCLUDE_DIR}")
include_directories ("${IMGUI_INCLUDE_DIR}")



#include_directories ("${CMAKE_SOURCE_DIR}/dependencies/glm/include")
#include_directories ("${CMAKE_SOURCE_DIR}/dependencies/glfw/include")

file(GLOB_RECURSE dep_sources ${SPIRVCROSS_INCLUDE_DIR}/*.cpp ${TINYOBJ_INCLUDE_DIR}/*.cc ${JSONCPP_INCLUDE_DIR}/*.cpp ${IMGUI_INCLUDE_DIR}/*.cpp)
file(GLOB_RECURSE legit_sources ./src/*.cpp)
file(GLOB_RECURSE legit_headers ./src/*.h)

add_executable(LegitEngine)

# Looks like without this the compiler still has too much optimization enabled in the debug build which makes it impossible to debug optimized code
target_compile_options(LegitEngine PRIVATE
  $<$<CONFIG:Debug>:-O0>
)

#source_group(TREE src/)
target_sources(LegitEngine PRIVATE ${dep_sources} ${legit_sources} ${legit_headers})
target_link_libraries(LegitEngine ${Vulkan_LIBRARIES} ${GLFW3_LIBRARY})
set_target_properties(LegitEngine PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}/bin/cmaked")
set_target_properties(LegitEngine PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/bin/cmake")