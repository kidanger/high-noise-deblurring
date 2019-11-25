#pragma once

//#include <algorithm>

#include "image.hpp"
#include "image_expr.hpp"

namespace utils {

    template <typename T>
    void transpose(img_t<T>& out, const img_t<T>& in)
    {
        if (&in == &out) {
            auto copy = in;
            return transpose(out, copy);
        }

        out.resize(in.h, in.w, in.d);
        for (int d = 0; d < in.d; d++) {
            for (int y = 0; y < in.h; y++) {
                for (int x = 0; x < in.w; x++) {
                    out(y, x, d) = in(x, y, d);
                }
            }
        }
    }

    template <typename T>
    img_t<T> add_padding(const img_t<T>& _f, int hw, int hh)
    {
        img_t<T> f(_f.w + hw*2, _f.h + hh*2, _f.d);
        f.set_value(T(0));
        slice(f, _(hw, -hw-1), _(hh, -hh-1)).map(_f);
        // replicate borders
        for (int y = 0; y < hh; y++) {
            for (int x = 0; x < f.w; x++) {
                for (int l = 0; l < f.d; l++) {
                    f(x, y, l) = f(x, 2*hh - y, l);
                    f(x, f.h-1-y, l) = f(x, f.h-1-2*hh+y, l);
                }
            }
        }
        for (int y = 0; y < f.h; y++) {
            for (int x = 0; x < hw; x++) {
                for (int l = 0; l < f.d; l++) {
                    f(x, y, l) = f(2*hw - x, y, l);
                    f(f.w-1-x, y, l) = f(f.w-1-2*hw+x, y, l);
                }
            }
        }
        return f;
    }

    template <typename T>
    img_t<T> add_padding(const img_t<T>& f, const img_t<T>& K)
    {
        return add_padding(f, K.w/2, K.h/2);
    }

    template <typename T>
    img_t<T> remove_padding(const img_t<T>& f, int hw, int hh)
    {
        return to_img(slice(f, _(hw, -hw-1), _(hh, -hh-1)));
    }

    template <typename T>
    img_t<T> remove_padding(const img_t<T>& f, const img_t<T>& K)
    {
        return remove_padding(f, K.w/2, K.h/2);
    }

}

