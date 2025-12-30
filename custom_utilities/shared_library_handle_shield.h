#if !defined(KRATOS_SHARED_LIBRARY_HANDLE_SHIELD_H_INCLUDED )
#define  KRATOS_SHARED_LIBRARY_HANDLE_SHIELD_H_INCLUDED

#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <mutex>

// --- Platform Abstraction Layer ---
#ifdef _WIN32
    #include <windows.h>
    #define PLATFORM_LOAD_LIB(name) (void*)LoadLibraryA(name)
    #define PLATFORM_FREE_LIB(handle) FreeLibrary((HMODULE)handle)
    #define PLATFORM_GET_SYM(handle, name) (void*)GetProcAddress((HMODULE)handle, name)
    #define PLATFORM_GET_ERROR() GetLastError()
#else
    #include <dlfcn.h>
    #define PLATFORM_LOAD_LIB(name) dlopen(name, RTLD_NOW | RTLD_GLOBAL)
    #define PLATFORM_FREE_LIB(handle) dlclose(handle)
    #define PLATFORM_GET_SYM(handle, name) dlsym(handle, name)
    #define PLATFORM_GET_ERROR() dlerror()
#endif

namespace Kratos
{

// forward declaration
namespace DLL
{
    inline std::string GetErrorMessage(
#ifdef _WIN32
        DWORD error
#else
        const char* error
#endif
    );
}

/**
 * Shield class to platform-independent load the shared library/DLL and unload when needed
 */
struct SharedLibraryHandleShield
{
    void* mHandle;
    std::string mlibraryName;

    SharedLibraryHandleShield(const std::string& LibraryName)
    {
        mlibraryName = LibraryName;
        mHandle = PLATFORM_LOAD_LIB(mlibraryName.c_str());

        if (!mHandle)
        {
#ifdef _WIN32
            DWORD error = PLATFORM_GET_ERROR();
            std::cerr << "CRITICAL ERROR: Could not load " << LibraryName
                      << ". Reason: " << DLL::GetErrorMessage(error) << std::endl;
#else
            const char* error = PLATFORM_GET_ERROR();
            std::cerr << "CRITICAL ERROR: Could not load " << LibraryName
                      << ". Reason: " << (error ? error : "Unknown error") << std::endl;
#endif
        }
    }

    ~SharedLibraryHandleShield()
    {
        if (mHandle)
        {
            PLATFORM_FREE_LIB(mHandle);
            std::cout << "Successfully unload shared library/DLL " << mlibraryName << std::endl;
        }
    }

    // Prevent accidental copying
    SharedLibraryHandleShield(const SharedLibraryHandleShield&) = delete;
    SharedLibraryHandleShield& operator=(const SharedLibraryHandleShield&) = delete;
};

namespace DLL
{

/**
 * Utility function to obtain the handle of the shared library using RAII approach
 */
inline void* GetSharedLibraryHandle(const std::string& LibraryName)
{
    // 1. Static map to keep track of all loaded libraries
    // 2. Mutex to make it thread-safe for multi-threaded calls
    static std::map<std::string, std::unique_ptr<SharedLibraryHandleShield>> s_Registry;
    static std::mutex s_Mutex;

    std::lock_guard<std::mutex> lock(s_Mutex);

    // Check if we already loaded this specific library
    auto it = s_Registry.find(LibraryName);
    if (it == s_Registry.end()) {
        // Not found: Create a new shield (this calls dlopen)
        s_Registry[LibraryName] = std::make_unique<SharedLibraryHandleShield>(LibraryName);
        return s_Registry[LibraryName]->mHandle;
    }

    // Found: Return the existing handle
    return it->second->mHandle;
}

/**
 * Get the specific symbol (function) within a handle (DLL)
 */
inline void* GetSymbol(void* handle, const std::string& SymbolName)
{
    return PLATFORM_GET_SYM(handle, SymbolName.c_str());
}

#ifdef _WIN32
inline DWORD GetError()
{
    return PLATFORM_GET_ERROR();
}
#else
inline const char* GetError()
{
    return PLATFORM_GET_ERROR();
}
#endif

inline std::string GetErrorMessage(
#ifdef _WIN32
    DWORD error
#else
    const char* error
#endif
)
{
#ifdef _WIN32
    if (error == 0) return "No error";

    LPSTR messageBuffer = nullptr;
    size_t size = FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL, error, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPSTR)&messageBuffer, 0, NULL);

    std::string message(messageBuffer, size);

    LocalFree(messageBuffer);

    return message + " (Code: " + std::to_string(error) + ")";
#else
    return (error) ? std::string(error) : "No error or error already cleared";
#endif
}

}

} // namespace Kratos

#undef PLATFORM_LOAD_LIB
#undef PLATFORM_FREE_LIB
#undef PLATFORM_GET_SYM
#undef PLATFORM_GET_ERROR

#endif // KRATOS_SHARED_LIBRARY_HANDLE_SHIELD_H_INCLUDED
