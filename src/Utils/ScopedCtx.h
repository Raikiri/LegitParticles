
template<typename T>
struct ScopedCtx
{
  ScopedCtx(entt::registry& reg, T&& val)
    : reg(reg)
  {
    reg.set<T>(std::forward<T>(val));
  }
  ScopedCtx(entt::registry& reg)
    : reg(reg)
  {
    reg.set<T>();
  }
  ~ScopedCtx()
  {
    reg.unset<T>();
  }
  entt::registry& reg;
};