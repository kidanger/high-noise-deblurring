#pragma once

#include <numeric>

#include "image.hpp"
#include "image_expr.hpp"
#include "vec2.hpp"
#include "smapa.h"

namespace optimization {

    template <typename T>
    struct Functional {
        typedef T value_type;
        auto gradient(const img_t<T>& x);
        auto prox(const img_t<T>& x, T tau);
        auto prox_adjoint(const img_t<T>& x, T tau);
    };

    template <typename T, typename T2>
    struct Operator {
        typedef T in_type;
        typedef T2 out_type;
        auto direct(const img_t<T>& x);
        auto adjoint(const img_t<T2>& x);
    };

    namespace operators {
        template <typename T>
        struct identity : Operator<T, T> {
            auto direct(const img_t<T>& x) {
                return to_expr(x);
            }
            auto adjoint(const img_t<T>& x) {
                return to_expr(x);
            }
        };

        template <typename T>
        struct gradient : Operator<T, vec2<T>> {
            img_t<vec2<T>> grad;
            img_t<T> div;
            gradient(const img_t<T> f) : grad(f.w, f.h, f.d), div(f.w, f.h, f.d) {}

            auto direct(const img_t<T>& x) {
                grad.gradients(x);
                return to_expr(grad);
            }

            auto adjoint(const img_t<vec2<T>>& x) {
                div.divergence(x);
                return -div;
            }
        };

        template <typename T>
        struct circular_gradient : Operator<T, vec2<T>> {
            img_t<vec2<T>> grad;
            img_t<T> div;
            circular_gradient(const img_t<T> f) : grad(f.w, f.h, f.d), div(f.w, f.h, f.d) {}

            auto direct(const img_t<T>& x) {
                grad.circular_gradients(x);
                return to_expr(grad);
            }

            auto adjoint(const img_t<vec2<T>>& x) {
                div.circular_divergence(x);
                return -div;
            }
        };

        template <typename T>
        struct circular_convolution : Operator<T, T> {
            img_t<std::complex<T>> ftK;
            img_t<std::complex<T>> ftx;

            circular_convolution(const img_t<T>& f, const img_t<T>& K) :
                    ftx(f.w, f.h, f.d),
                    ftK(f.w, f.h, f.d) {
                assert(K.w % 2);
                assert(K.h % 2);

                ftK.padcirc(K);
                ftK.map(ftK * T(K.d) / K.sum());
                ftK.fft(ftK);
            }

            auto direct(const img_t<T>& x) {
                ftx.map(x);
                ftx.fft(ftx);
                ftx.map(ftx * ftK);
                ftx.ifft(ftx);
                return std::real(ftx);
            }

            auto adjoint(const img_t<T>& x) {
                ftx.map(x);
                ftx.fft(ftx);
                ftx.map(ftx * std::conj(ftK));
                ftx.ifft(ftx);
                return std::real(ftx);
            }
        };
    }

    namespace functionals {
        template <typename T>
        struct indicator_always_zero : Functional<T> {
            auto prox(const img_t<T>& x, T tau) {
                return to_expr(x);
            }
        };

        template <typename T>
        struct xminusyL2 : Functional<T> {
            img_t<T> y;
            xminusyL2(const img_t<T>& y) : y(y) {}

            auto gradient(const img_t<T>& x) {
                return x - y;
            }
        };

        template <typename T>
        struct KxminusyL2 : Functional<T> {
            img_t<T> y;
            img_t<T> diff;
            optimization::operators::circular_convolution<T> conv;
            KxminusyL2(const img_t<T>& y, const img_t<T>& K) : y(y), conv(y, K), diff(y) {}

            auto gradient(const img_t<T>& x) {
                diff.map(conv.direct(x) - y);
                return conv.adjoint(diff);
            }
        };

        template <typename T>
        struct lambda_norm1_2 : Functional<vec2<T>> {
            T lambda;
            lambda_norm1_2(T lambda=T(1)) : lambda(lambda) {}

            auto prox_adjoint(const img_t<vec2<T>>& x, T tau) {
                return x / std::max(std::hypot(x) / lambda, T(1));
            }
        };

        template <typename T>
        struct nlgradients_norm2 : Functional<T> {
            struct point_t {
                T w;
                int x;
                int y;
                int d;
            };
            img_t<std::vector<point_t>> weights;
            img_t<T> nlu;
            int S, N;
            T h;

            nlgradients_norm2(const img_t<T>& u, int S, int N, T h) : nlu(u), weights(u.w, u.h, u.d), S(S), N(N), h(h) {
                reweight(u);
            }

            void reweight(const img_t<T>& u) {
                int hS = S / 2;
                int hN = N / 2;
                T h2 = std::pow(h, T(2));
                for (int d = 0; d < u.d; d++) {
                    printf("d=%d\n", d);
                    for (int y = 0; y < u.h; y++) {
                        for (int x = 0; x < u.w; x++) {
                            weights(x, y, d).resize(0);
                            // search window
                            for (int pdy = -hS; pdy <= hS; pdy++) {
                                for (int pdx = -hS; pdx <= hS; pdx++) {
                                    if (!u.inside(x + pdx, y + pdy))
                                        continue;
                                    // window similarity
                                    T similarity(0);
                                    for (int wdy = -hN; wdy <= hN; wdy++) {
                                        for (int wdx = -hN; wdx <= hN; wdx++) {
                                            if (u.inside(x + wdx, y + wdy) && u.inside(x + pdx + wdx, y + pdy + wdy))
                                                similarity += std::pow(u(x + wdx, y + wdy, d) - u(x + pdx + wdx, y + pdy + wdy, d), T(2));
                                        }
                                    }
                                    //similarity -= 2*variance;
                                    similarity = std::exp(- similarity / h2);
                                    weights(x, y, d).push_back(point_t{.w=similarity, .x=x+pdx, .y=y+pdy, .d=d});
                                }
                            }
                        }
                    }
                }

                for (int i = 0; i < weights.size; i++) {
                    T sum = 0;
                    for (auto& p : weights[i])
                        sum += p.w;
                    for (auto& p : weights[i])
                        p.w /= sum;
                }
            }

            auto nl(const img_t<T>& u) {
                nlu.similar(u);
                nlu.set_value(T(0));
                for (int i = 0; i < nlu.size; i++) {
                    for (auto& p : weights[i])
                        nlu[i] += p.w * u(p.x, p.y, p.d);
                }
                return nlu;
            }

            auto gradient(const img_t<T>& u) {
                return u - nl(u);
            }
        };
    }

    namespace methods {
        SMART_PARAMETER_STR(DEBUG_OUT);

        template <typename T, class P>
        void gradient_descent(img_t<T>& u, P& p, float step, int iter) {
            for (int i = 0; i < iter; i++) {
                u.map(u - step * p.gradient(u));
                fprintf(stderr, "step %d/%d\n", i+1, iter);
                if (DEBUG_OUT())
                    u.save(DEBUG_OUT());
            }
        }

        template <typename T, class P>
        void condat(img_t<T>& x, P& p, int iter, T rho, T tau, T sigma) {
            img_t<T> xn(x);
            img_t<T> xh(x);
            img_t<vec2<T>> u(x.w, x.h, x.d);
            img_t<vec2<T>> uh(x.w, x.h, x.d);
            img_t<T> g(x);

            u.set_value(vec2<T>(0, 0));
            for (int i = 0; i < iter; i++) {
                xh.map(x - tau * (p.f.gradient(x) + p.L.adjoint(u)));

                xh.map(p.g.prox(xh, tau));

                xn.map(rho * xh + (1.f - rho) * x);

                g.map(T(2) * xh - x);
                uh.map(u + sigma * p.L.direct(g));

                uh.map(p.h.prox_adjoint(uh, sigma));

                u.map(rho * uh + (1.f - rho) * u);

                x.copy(xn);
            }
        }
    }
}

