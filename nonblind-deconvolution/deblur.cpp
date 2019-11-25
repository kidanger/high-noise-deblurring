#include <iostream>
#include <random>

#include "smapa.h"
#include "deblur.hpp"
#include "utils.hpp"
#include "edgetaper.hpp"

SMART_PARAMETER(int, NORMALIZE_INPUT, 1);
SMART_PARAMETER(int, EDGETAPER, 1);

SMART_PARAMETER(float, CONTINUATION_BETA_INIT, 1);
SMART_PARAMETER(float, CONTINUATION_BETA_RATE, 2.f * std::sqrt(2.f));
SMART_PARAMETER(float, CONTINUATION_BETA_MAX, std::pow(2.f, 8.f));

int main(int argc, char** argv) {
    if (argc != 5)
        return std::cerr << "Usage: " << argv[0] << " input kernel output lambda" << std::endl, 1;

    img_t<float> f = img_t<float>::load(argv[1]);
    img_t<float> K = img_t<float>::load(argv[2]);
    K.map(K / K.sum());
    img_t<float> u;

    float max = f.max();
    if (NORMALIZE_INPUT())
        f.map(f / f.max());

    if (EDGETAPER()) {
        f = utils::add_padding(f, K);
        edgetaper(f, f, K, 3);
    }

    float lambda = atof(argv[4]);
    deblur::rof::split_continuation(u, f, K, 2.f / lambda, CONTINUATION_BETA_INIT(),
                                    CONTINUATION_BETA_RATE(), CONTINUATION_BETA_MAX());

    if (EDGETAPER()) {
        u = utils::remove_padding(u, K);
    }

    if (NORMALIZE_INPUT())
        u.map(u * max);

    u.save(argv[3]);
}

