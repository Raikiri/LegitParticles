#include <string>
#include <fstream>
#include <functional>

namespace legit
{
  class Core;
  class ShaderModule
  {
  public:
    vk::ShaderModule GetHandle()
    {
      return shaderModule.get();
    }
    uint32_t GetHash()
    {
      return hash;
    }
    ShaderModule(vk::Device device, const std::vector<uint32_t> &bytecode)
    {
      Init(device, bytecode);
    }
  private:
    void Init(vk::Device device, const std::vector<uint32_t> &bytecode)
    {
      auto shaderModuleCreateInfo = vk::ShaderModuleCreateInfo()
        .setCodeSize(bytecode.size() * sizeof(uint32_t))
        .setPCode(bytecode.data());
      this->shaderModule = device.createShaderModuleUnique(shaderModuleCreateInfo);
      this->hash = 0;
      for(auto b : bytecode)
      {
        this->hash ^= b; //actually terrible hash
      }
    }
    vk::UniqueShaderModule shaderModule;
    uint32_t hash;
    friend class Core;
  };
}