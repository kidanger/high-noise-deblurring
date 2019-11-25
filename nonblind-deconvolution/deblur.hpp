#pragma once

#include "image.hpp"
#include "vec2.hpp"
#include "optimization.hpp"
#include "fft.hpp"
#include "utils.hpp"

namespace deblur {
    namespace rof {

        template <typename T>
        void split_continuation(img_t<T>& x, const img_t<T>& f, const img_t<T>& K,
                                T lambda, T beta_init, T beta_rate, T beta_max) {
            assert(K.w % 2);
            assert(K.h % 2);
            optimization::operators::gradient<T> gradient(f);
            using gradtype = typename decltype(gradient)::out_type;

            img_t<gradtype> w(f.w, f.h, f.d);

            img_t<std::complex<T>> w1_ft(f.w, f.h, f.d);
            img_t<std::complex<T>> x_ft(f.w, f.h, f.d);
            img_t<std::complex<T>> f_ft(f.w, f.h, f.d);
            f_ft.map(f);
            f_ft.fft(f_ft);

            img_t<std::complex<T>> dx_otf(f.w, f.h, f.d);
            img_t<T> dx(3, 3, f.d);
            for (int l = 0; l < f.d; l++) {
                dx(0, 1, l) = 0; dx(1, 1, l) = -1; dx(2, 1, l) = 1;
            }
            dx_otf.padcirc(dx);
            dx_otf.fft(dx_otf);

            img_t<std::complex<T>> dy_otf(f.w, f.h, f.d);
            img_t<T> dy(3, 3, f.d);
            for (int l = 0; l < f.d; l++) {
                dy(1, 0, l) = 0; dy(1, 1, l) = -1; dy(1, 2, l) = 1;
            }
            dy_otf.padcirc(dy);
            dy_otf.fft(dy_otf);

            img_t<std::complex<T>> K_otf(f.w, f.h, f.d);
            K_otf.padcirc(K);
            K_otf.map(K_otf * std::complex<T>(K.d) / K.sum());
            K_otf.fft(K_otf);

            auto Ktf = std::conj(K_otf) * f_ft;

            img_t<T> KtK(f.w, f.h, f.d);
            KtK.map(std::pow(std::abs(K_otf), T(2)));

            img_t<T> DtD(f.w, f.h, f.d);
            auto Fdx = std::pow(std::abs(dx_otf), T(2));
            auto Fdy = std::pow(std::abs(dy_otf), T(2));
            DtD.map(Fdx + Fdy);

            T beta = beta_init;
            x = f;
            while (beta < beta_max) {
                T gamma = beta / lambda;
                auto denom = KtK + gamma*DtD;

                for (int inner = 0; inner < 1; inner++) {
                    auto gx = gradient.direct(x);
                    //w.map(std::max(std::abs(gx) - T(1) / beta, T(0)) * std::sign(gx));  // anisotropic
                    w.map(gx * (T(1) - T(1) / (std::max(T(1), beta * std::hypot(gx)))));  // isotropic

                    w1_ft.map(gradient.adjoint(w));
                    w1_ft.fft(w1_ft);
                    auto num = Ktf + gamma * w1_ft;
                    x_ft.map(num / denom);
                    x_ft.ifft(x_ft);
                    x.map(std::real(x_ft));
                }

                beta *= beta_rate;
            }
        }
    }
};

