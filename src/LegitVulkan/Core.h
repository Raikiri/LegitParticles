#include <iostream>

#include <vector>
#include <set>
namespace legit
{
  class Swapchain;
  class RenderGraph;
  class Core
  {
  public:    
    inline Core(const char** instanceExtensions, uint32_t instanceExtensionsCount, WindowDesc* compatibleWindowDesc, bool enableDebugging);

    inline ~Core();
    inline void ClearCaches();
    inline std::unique_ptr<Swapchain> CreateSwapchain(WindowDesc windowDesc, glm::uvec2 defaultSize, uint32_t imagesCount, vk::PresentModeKHR preferredMode);


    template<typename Handle, typename Loader>
    static inline void SetObjectDebugName(vk::Device logicalDevice, Loader& loader, Handle objHandle, std::string name)
    {
      auto nameInfo = vk::DebugUtilsObjectNameInfoEXT()
        .setObjectHandle(uint64_t(typename Handle::CType(objHandle)))
        .setObjectType(objHandle.objectType)
        .setPObjectName(name.c_str());
      if(loader.vkSetDebugUtilsObjectNameEXT)
        auto res = logicalDevice.setDebugUtilsObjectNameEXT(&nameInfo, loader);
    }

    template<typename T>
    void SetObjectDebugName(T objHandle, std::string name)
    {
      SetObjectDebugName(logicalDevice.get(), loader, objHandle, name);
    }
    
    inline void SetDebugName(legit::ImageData* imageData, std::string name);

    inline vk::CommandPool GetCommandPool();
    inline vk::Device GetLogicalDevice();
    inline vk::PhysicalDevice GetPhysicalDevice();
    inline legit::RenderGraph* GetRenderGraph();
    inline std::vector<vk::UniqueCommandBuffer> AllocateCommandBuffers(size_t count);
    inline vk::UniqueSemaphore CreateVulkanSemaphore();
    inline vk::UniqueFence CreateFence(bool state);
    inline void WaitForFence(vk::Fence fence);
    inline void ResetFence(vk::Fence fence);
    inline void WaitIdle();
    inline vk::Queue GetGraphicsQueue();
    inline vk::Queue GetPresentQueue();
    inline uint32_t GetDynamicMemoryAlignment();
    inline legit::DescriptorSetCache* GetDescriptorSetCache();
    inline legit::PipelineCache* GetPipelineCache();
  private:

    static inline vk::UniqueInstance CreateInstance(const std::vector<const char*>& instanceExtensions, const std::vector<const char*>& validationLayers);

    static inline vk::UniqueHandle<vk::DebugUtilsMessengerEXT, vk::detail::DispatchLoaderDynamic> CreateDebugUtilsMessenger(vk::Instance instance, vk::PFN_DebugUtilsMessengerCallbackEXT debugCallback, vk::detail::DispatchLoaderDynamic& loader);

    static inline vk::Bool32 DebugMessageCallback(
      vk::DebugUtilsMessageSeverityFlagBitsEXT       messageSeverity,
      vk::DebugUtilsMessageTypeFlagsEXT              messageTypes,
      const vk::DebugUtilsMessengerCallbackDataEXT * pCallbackData,
      void *pUserData );

    friend class Swapchain;
    static inline vk::PhysicalDevice FindPhysicalDevice(vk::Instance instance);

    static inline QueueFamilyIndices FindQueueFamilyIndices(vk::PhysicalDevice physicalDevice, vk::SurfaceKHR surface);

    static inline vk::UniqueDevice CreateLogicalDevice(vk::PhysicalDevice physicalDevice, QueueFamilyIndices familyIndices, std::vector<const char*> deviceExtensions, std::vector<const char*> validationLayers);
    static inline vk::Queue GetDeviceQueue(vk::Device logicalDevice, uint32_t queueFamilyIndex);
    static inline vk::UniqueCommandPool CreateCommandPool(vk::Device logicalDevice, uint32_t familyIndex);

    vk::UniqueInstance instance;
    vk::detail::DispatchLoaderDynamic loader;
    vk::UniqueHandle<vk::DebugUtilsMessengerEXT, vk::detail::DispatchLoaderDynamic> debugUtilsMessenger;
    vk::PhysicalDevice physicalDevice;
    vk::UniqueDevice logicalDevice;
    vk::UniqueCommandPool commandPool;
    vk::Queue graphicsQueue;
    vk::Queue presentQueue;

    std::unique_ptr<legit::DescriptorSetCache> descriptorSetCache;
    std::unique_ptr<legit::PipelineCache> pipelineCache;
    std::unique_ptr<legit::RenderGraph> renderGraph;


    QueueFamilyIndices queueFamilyIndices;
  };
}