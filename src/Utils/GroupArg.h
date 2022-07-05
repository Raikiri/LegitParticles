namespace entt
{
  template<typename ...Components>
  using group_type = entt::group<entt::exclude_t<>, entt::get_t<>, Components...>;
}