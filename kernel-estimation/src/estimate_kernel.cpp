#include <iostream>
#include <random>

#include "estimate_kernel.hpp"

#include "args.hxx"
#include <iostream>

static options parse_args(int argc, char** argv)
{
    args::ArgumentParser parser("");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<int> ks(parser, "ks", "blur kernel size", args::Options::Required);
    args::Positional<std::string> input(parser, "input", "input blurry image file", args::Options::Required);
    args::Positional<std::string> output(parser, "output", "kernel output file", args::Options::Required);
    args::ValueFlag<float> lambda(parser, "lambda", "L0 regularization weight", {"lambda"}, 4e-3f);
    args::ValueFlag<float> lambda_ratio(parser, "lambda-ratio", "decay of lambda", {"lambda-ratio"}, 1/1.1f);
    args::ValueFlag<float> lambda_min(parser, "lambda-min", "L0 regularization weight minimum value", {"lambda-min"}, 1e-2f);
    args::ValueFlag<float> gamma(parser, "gamma", "kernel regularization weight", {"gamma"}, 20.f);
    args::ValueFlag<int> iterations(parser, "iterations", "number of iterations per scale", {"iterations"}, 2);
    args::Flag multiscale(parser, "no-multiscale", "disable the multiscale scheme", {"no-multiscale"});
    args::ValueFlag<float> scalefactor(parser, "scale-factor", "downsampling factor", {"scale-factor"}, 0.5f);
    args::ValueFlag<float> kernel_threshold_max(parser, "kernel-threshold-max", "threshold the kernel at max(kernel)*kernel-threshold-max ", {"kernel-threshold-max"}, 0.f);
    args::ValueFlag<bool> remove_isolated(parser, "remove-isolated", "remove isolated connected component of the kernel", {"remove-isolated"}, false);
    args::ValueFlag<std::string> outputsharp(parser, "output-sharp", "output the sharp image to file", {"output-sharp"});
    args::ValueFlag<std::string> debug(parser, "debug", "output all kernels, sharp and blurry images", {"debug"});
    args::Flag better_kernel(parser, "no-better-kernel-estimation", "...", {"better-kernel-estimation"}, false);
    args::Flag warmg(parser, "no-warmg", "...", {"warmg"}, false);
    args::Flag warmk(parser, "warmk", "...", {"warmk"}, false);
    args::ValueFlag<float> upscaleblur(parser, "upscale-blur", "", {"upscale-blur"}, 0);
    args::ValueFlag<float> downscaleblur(parser, "downscale-blur", "", {"downscale-blur"}, 1.6f);
    args::Flag verbose(parser, "verbose", "output more information", {"verbose"});
    args::ValueFlag<std::string> initu(parser, "initu", "", {"initu"});
    args::Flag admmu(parser, "admmu", "", {"admmu"});
    args::ValueFlag<float> admm_mu(parser, "admm-mu", "", {"admm-mu"});
    args::ValueFlag<float> k_l1(parser, "k-l1", "", {"k-l1"}, 0.5f);
    args::ValueFlag<bool> use_filters(parser, "use-filters", "", {"use-filters"}, false);

    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        exit(0);
    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    } catch (const args::ValidationError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    }

    options opts;
    opts.ks = args::get(ks);
    opts.input = args::get(input);
    opts.output = args::get(output);
    opts.lambda = args::get(lambda);
    opts.lambda_ratio = args::get(lambda_ratio);
    opts.lambda_min = args::get(lambda_min);
    opts.gamma = args::get(gamma);
    opts.iterations = args::get(iterations);
    opts.multiscale = !args::get(multiscale);
    opts.scalefactor = args::get(scalefactor);
    opts.kernel_threshold_max = args::get(kernel_threshold_max);
    opts.remove_isolated = args::get(remove_isolated);
    opts.outputsharp = args::get(outputsharp);
    opts.verbose = args::get(verbose);
    opts.debug = args::get(debug);
    opts.better_kernel = !args::get(better_kernel);
    opts.warmg = !args::get(warmg);
    opts.warmk = args::get(warmk);
    opts.upscaleblur = args::get(upscaleblur);
    opts.downscaleblur = args::get(downscaleblur);
    opts.initu = args::get(initu);
    opts.admmu = args::get(admmu);
    opts.admmu_mu = args::get(admm_mu);
    opts.k_l1 = args::get(k_l1);
    opts.use_filters = args::get(use_filters);

    if (getenv("PRINT_PARAMS") && *getenv("PRINT_PARAMS") == '1') {
        printf("params\t");
        printf("%d", opts.ks);
        printf(",%f", opts.lambda);
        printf(",%f", opts.lambda_ratio);
        printf(",%f", opts.lambda_min);
        printf(",%f", opts.gamma);
        printf(",%d", opts.iterations);
        printf(",%d", opts.multiscale);
        printf(",%f", opts.scalefactor);
        printf(",%f", opts.kernel_threshold_max);
        printf(",%d", opts.remove_isolated);
        printf(",%d", opts.better_kernel);
        printf(",%d", opts.warmg);
        printf(",%d", opts.warmk);
        printf(",%f", opts.upscaleblur);
        printf(",%f", opts.downscaleblur);
    }

    return opts;
}

int main(int argc, char** argv) {
    struct options opts = parse_args(argc, argv);

    img_t<float> v = img_t<float>::load(opts.input);

    preprocess_image(v, v, opts);

    img_t<float> initu;
    if (opts.initu.empty()) {
        initu = v;
    } else {
        initu = img_t<float>::load(opts.initu);
        preprocess_image(initu, initu, opts);
        if (opts.multiscale) {
            fprintf(stderr, "wait, both initu and multiscale are enabled\n");
        }
    }

    img_t<float> k;
    img_t<float> u;
    if (opts.multiscale) {
        multiscale_l0_kernel_estimation(k, u, v, opts);
    } else {
        l0_kernel_estimation(k, u, v, initu, opts);
    }

    k.save(opts.output);
    if (!opts.outputsharp.empty()) {
        u.save(opts.outputsharp);
    }
}

