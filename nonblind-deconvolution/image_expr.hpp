#pragma once

#include <type_traits>

template <typename T>
class img_t;

class dummy_img_expr_t {
};

template <typename T>
class img_expr_t : public dummy_img_expr_t {
public:
    typedef T value_type;
    //int w, h, d;

    // not defined here
    // T operator[](int i) const;
};

template <typename T>
class scalar_img_expr_t : public img_expr_t<T> {
public:
    T val;
    int size = 0;
    int w = 0;
    int h = 0;
    int d = 0;

    scalar_img_expr_t(T val) : val(val) {
    }

    T operator[](int i) const {
        return val;
    }

    template <typename E>
    bool similar(const E& o) const {
        return true;
    }
};

template <typename T>
class raw_img_expr_t : public img_expr_t<T> {
public:
    const img_t<T>* img;
    int size, w, h, d;

    raw_img_expr_t(const img_t<T>* img) : img(img), size(img->size),
                                          w(img->w), h(img->h), d(img->d) {
    }

    T operator[](int i) const {
        return (*img)[i];
    }

    template <typename E>
    bool similar(const E& o) const {
        return o.similar(*img);
    }
};

template <typename E>
typename std::enable_if<std::is_base_of<dummy_img_expr_t, E>::value, E>::type
to_expr(const E& e) {
    return e;
}

// important to avoid unecessary image copy
// use code should use to_expr when returning a raw image
template <typename T>
raw_img_expr_t<T> to_expr(const img_t<T>& i) {
    return raw_img_expr_t<T>(&i);
}

template <typename T>
typename std::enable_if<!std::is_base_of<dummy_img_expr_t, T>::value, scalar_img_expr_t<T>>::type
to_expr(T t) {
    return scalar_img_expr_t<T>(t);
}


template <typename T>
class func0_img_expr_t : public img_expr_t<T> {
public:
    std::function<T(int)> f;
    int size = 0;
    int w = 0;
    int h = 0;
    int d = 0;

    func0_img_expr_t(std::function<T(int i)> f) : f(f) {
    }

    T operator[](int i) const {
        return f(i);
    }

    template <typename T2>
    bool similar(const img_t<T2>& o) const {
        return true;
    }
};

template <typename T, typename E>
class func1_img_expr_t : public img_expr_t<T> {
public:
    std::function<T(typename E::value_type)> f;
    E e;
    int size, w, h, d;

    func1_img_expr_t(std::function<T(typename E::value_type)> f, const E& e)
        : f(f), e(e), size(e.size), w(e.w), h(e.h), d(e.d) {
    }

    T operator[](int i) const {
        return f(e[i]);
    }

    template <typename E2>
    bool similar(const E2& o) const {
        return e.similar(o);
    }
};

template <typename T, typename E1, typename E2>
class func2_img_expr_t : public img_expr_t<T> {
public:
    std::function<T(typename E1::value_type, typename E2::value_type)> f;
    E1 e1;
    E2 e2;
    int size, w, h, d;

    func2_img_expr_t(std::function<T(typename E1::value_type, typename E2::value_type)> f,
                     const E1& e1, const E2& e2) : f(f), e1(e1), e2(e2), size(e1.size | e2.size),
                                                   w(e1.w | e2.w), h(e1.h | e2.h), d(e1.d | e2.d) {
    }

    T operator[](int i) const {
        return f(e1[i], e2[i]);
    }

    template <typename E3>
    bool similar(const E3& o) const {
        return e1.similar(o) && e2.similar(o);
    }
};

#define DEFINE_EXPR_1(name, operand) \
template <typename E> \
class name : \
    public img_expr_t<decltype(operand std::declval<typename E::value_type>())> { \
    typedef decltype(operand std::declval<typename E::value_type>()) res_type; \
public: \
    E e; \
    int size, w, h, d; \
\
    name(const E& e) : e(e), size(e.size), w(e.w), h(e.h), d(e.d) { \
    } \
\
    res_type operator[](int i) const { \
        return operand e[i]; \
    } \
    \
    template <typename E3> \
    bool similar(const E3& o) const { \
        return e.similar(o); \
    } \
}; \
\
template <typename E> \
name<decltype(to_expr(std::declval<E>()))> operator operand(const E& e) { \
    return name<decltype(to_expr(e))>(to_expr(e)); \
}; \

#define DEFINE_EXPR_2(name, operand) \
template <typename E1, typename E2> \
class name : \
    public img_expr_t<decltype(std::declval<typename E1::value_type>() operand std::declval<typename E2::value_type>())> { \
    typedef decltype(std::declval<typename E1::value_type>() operand std::declval<typename E2::value_type>()) res_type;\
public: \
    E1 e1; \
    E2 e2; \
    int size, w, h, d; \
\
    name(const E1& e1, const E2& e2) : e1(e1), e2(e2), size(e1.size | e2.size), \
                                       w(e1.w | e2.w), h(e1.h | e2.h), d(e1.d | e2.d) { \
    } \
\
    res_type operator[](int i) const { \
        return e1[i] operand e2[i]; \
    } \
    \
    template <typename E3> \
    bool similar(const E3& o) const { \
        return e1.similar(o) && e2.similar(o); \
    } \
}; \
\
template <typename E1, typename E2> \
name<decltype(to_expr(std::declval<E1>())), decltype(to_expr(std::declval<E2>()))> operator operand(const E1& e1, const E2& e2) { \
    return name<decltype(to_expr(e1)), decltype(to_expr(e2))>(to_expr(e1), to_expr(e2)); \
}; \

DEFINE_EXPR_1(minus_img_expr_t, -)
DEFINE_EXPR_2(add_img_expr_t, +)
DEFINE_EXPR_2(sub_img_expr_t, -)
DEFINE_EXPR_2(div_img_expr_t, /)
DEFINE_EXPR_2(mul_img_expr_t, *)

#define DEFINE_FUNC_1(name, call) \
    template <typename E> \
    auto name(const E& e) { \
        return func1_img_expr_t<decltype(call(std::declval<typename decltype(to_expr(e))::value_type>())), \
                                decltype(to_expr(e))>( \
           [](auto e) { return call(e); }, to_expr(e)); \
    } \

#define DEFINE_FUNC_2(name, call) \
    template <typename E1, typename E2> \
    auto name(const E1& e1, const E2& e2) { \
        return func2_img_expr_t<decltype(call(std::declval<typename decltype(to_expr(e1))::value_type>(), \
                                              std::declval<typename decltype(to_expr(e2))::value_type>())), \
                                decltype(to_expr(e1)), decltype(to_expr(e2))> \
            ([](auto e1, auto e2) { return call(e1, e2); }, to_expr(e1), to_expr(e2)); \
    } \

#include "better_than_std.hpp"
#include "vec2.hpp"

namespace std {
    DEFINE_FUNC_1(conj, std::conj)
    DEFINE_FUNC_1(real, std::real)
    DEFINE_FUNC_1(hypot, std::hypot)
    DEFINE_FUNC_1(log, std::log)
    DEFINE_FUNC_1(abs, std::abs)
    DEFINE_FUNC_1(arg, std::arg)
    DEFINE_FUNC_1(exp, std::exp)
    DEFINE_FUNC_1(floor, std::floor)
    DEFINE_FUNC_1(round, std::round)

    DEFINE_FUNC_2(max, std::max_noref)
    DEFINE_FUNC_2(min, std::min_noref)
    DEFINE_FUNC_2(hypot, std::hypot)
    DEFINE_FUNC_2(pow, std::pow)

    DEFINE_FUNC_1(sign, std::sgn)
}

template <typename T, typename E>
img_t<T> to_img(const E& e)
{
    img_t<T> img(e.w, e.h, e.d);
    img.map(e);
    return img;
}
template <typename T>
const img_t<T>& to_img(const img_t<T>& img)
{
    return img;
}
template <typename T, typename T2>
const img_t<T2>& to_img(const img_t<T>& img)
{
    img_t<T2> img2; img2.resize(img); img2.map(img);
    return img2;
}
template <typename E>
img_t<typename E::value_type> to_img(const E& e)
{
    return to_img<typename E::value_type>(e);
}

static struct slice_t {
    int s, e;
    slice_t operator()(int s) const {
        return {.s=s, .e=s};
    }
    slice_t operator()(int s, int e) const {
        return {.s=s, .e=e};
    }
} _ = {.s=0, .e=-1};

template <typename E>
class slice_img_expr_t : public img_expr_t<typename E::value_type> {
    using T = typename E::value_type;
public:
    E e;
    slice_t _w, _h, _d;
    int size, w, h, d;

    slice_img_expr_t(const E& e, slice_t _w, slice_t _h, slice_t _d)
        : e(e), _w(_w), _h(_h), _d(_d) {
        this->_w.s = (_w.s + e.w) % e.w;
        this->_h.s = (_h.s + e.h) % e.h;
        this->_d.s = (_d.s + e.d) % e.d;
        this->_w.e = (_w.e + e.w) % e.w;
        this->_h.e = (_h.e + e.h) % e.h;
        this->_d.e = (_d.e + e.d) % e.d;
        w = this->_w.e - this->_w.s + 1;
        h = this->_h.e - this->_h.s + 1;
        d = this->_d.e - this->_d.s + 1;
        size = w * h * d;
        //assert(similar(e));
        // check other size compatibility with e
    }

    T operator[](int i) const {
        int dd = i % d;
        int xy = i / d;
        int x = xy % w;
        int y = xy / w;
        x += _w.s;
        y += _h.s;
        dd += _d.s;
        int j = dd + e.d * (x + e.w * y);
        return e[j];
    }

    template <typename E2>
    bool similar(const E2& o) const {
        return o.similar(*this);
    }
};

template <typename T>
class mappable_slice_img_expr_t : public img_expr_t<T> {
public:
    img_t<T>* img;
    slice_t _w, _h, _d;
    int size, w, h, d;

    mappable_slice_img_expr_t(img_t<T>* img, slice_t _w, slice_t _h, slice_t _d)
        : img(img), _w(_w), _h(_h), _d(_d) {
        this->_w.s = (_w.s + img->w) % img->w;
        this->_h.s = (_h.s + img->h) % img->h;
        this->_d.s = (_d.s + img->d) % img->d;
        this->_w.e = (_w.e + img->w) % img->w;
        this->_h.e = (_h.e + img->h) % img->h;
        this->_d.e = (_d.e + img->d) % img->d;
        w = this->_w.e - this->_w.s + 1;
        h = this->_h.e - this->_h.s + 1;
        d = this->_d.e - this->_d.s + 1;
        size = w * h * d;
        // check other size compatibility with img
    }

    T operator[](int i) const {
        int dd = i % d;
        int xy = i / d;
        int x = xy % w;
        int y = xy / w;
        x += _w.s;
        y += _h.s;
        dd += _d.s;
        int j = dd + img->d * (x + img->w * y);
        return (*img)[j];
    }

    T& operator[](int i) {
        int dd = i % d;
        int xy = i / d;
        int x = xy % w;
        int y = xy / w;
        x += _w.s;
        y += _h.s;
        dd += _d.s;
        int j = dd + img->d * (x + img->w * y);
        return (*img)[j];
    }

    template <typename E>
    void map(const E& o) {
        auto e = to_expr(o);
        assert(e.similar(*this));
        for (int i = 0; i < size; i++)
            (*this)[i] = e[i];
    }

    template <typename E>
    bool similar(const E& o) const {
        bool similar = w == o.w && h == o.h && d == o.d;
        if (!similar)
            fprintf(stderr, "%dx%dx%d (type %s) != %dx%dx%d (type %s)\n",
                    w, h, d, typeid(*this).name(), o.w, o.h, o.d, typeid(o).name());
        return similar;
    }
};

template <typename E>
auto slice(const E& e, slice_t w, slice_t h, slice_t d=_)
{
    return slice_img_expr_t<decltype(to_expr(e))>(to_expr(e), w, h, d);
}

template <typename T>
auto slice(img_t<T>& i, slice_t w, slice_t h, slice_t d=_)
{
    return mappable_slice_img_expr_t<T>(&i, w, h, d);
}

