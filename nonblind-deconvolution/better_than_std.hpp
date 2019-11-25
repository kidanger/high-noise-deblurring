#pragma once

namespace std {
    template <typename T> constexpr T sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
}

namespace std {
    // TODO: find a way to fix that

    template<typename T> constexpr
    T min_noref(const T& a, const T& b) {
        return a > b ? b : a;
    }

    template<typename T> constexpr
    T max_noref(const T& a, const T& b) {
        return a > b ? a : b;
    }
}

