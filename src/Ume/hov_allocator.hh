#ifndef HOV_ALLOCATOR_HH
#define HOV_ALLOCATOR_HH

#include <cstddef>
#include <new>

extern "C" {
#include "hov.h"
#include "hov_alloc.h"
}

template <class T>
struct HovAllocator {
    typedef T value_type;

    HovAllocator() = default;
    
    template <class U> 
    constexpr HovAllocator(const HovAllocator<U>&) noexcept {}

    T* allocate(std::size_t n) {
        if (n > std::size_t(-1) / sizeof(T)) {
            throw std::bad_alloc();
        }
        if (auto p = static_cast<T*>(hov_alloc_data(n * sizeof(T)))) {
            return p;
        }
        throw std::bad_alloc();
    }

    void deallocate(T* p, std::size_t n) noexcept {
        (void)n;
        hov_free_data(p);
    }
};

template <class T, class U>
bool operator==(const HovAllocator<T>&, const HovAllocator<U>&) { return true; }

template <class T, class U>
bool operator!=(const HovAllocator<T>&, const HovAllocator<U>&) { return false; }

#endif // HOV_ALLOCATOR_HH
