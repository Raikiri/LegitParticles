namespace legit
{
    struct WindowDesc
    {
      #if defined(WIN32)
        HINSTANCE hInstance;
        HWND hWnd;
      #else
        wl_display *display;
        wl_surface *surface;
      #endif
    };

    #if defined(WIN32)
      static vk::UniqueSurfaceKHR CreateSurface(vk::Instance instance, WindowDesc desc)
      {
        auto surfaceCreateInfo = vk::Win32SurfaceCreateInfoKHR()
          .setHwnd(desc.hWnd)
          .setHinstance(desc.hInstance);

        return instance.createWin32SurfaceKHRUnique(surfaceCreateInfo);

      }
    #else
      static vk::UniqueSurfaceKHR CreateSurface(vk::Instance instance, WindowDesc desc)
      {
        auto surfaceCreateInfo = vk::WaylandSurfaceCreateInfoKHR()
          .setDisplay(desc.display)
          .setSurface(desc.surface);

        return instance.createWaylandSurfaceKHRUnique(surfaceCreateInfo);

      }
    #endif
}