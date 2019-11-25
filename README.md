
Usage
=====

Use the `run.sh` bash script to run the kernel estimation, denoising and deblurring.

At the strict minimum, you can give it your image as the first argument and it will estimate the noise level and produce the result alongside of your image:
```bash
$ bash run.sh example.png
# this produces the file "example.png-deblurred.png"
```

Other arguments (in order):

2. kernel size. Should be an odd number. 35 by default.
3. noise standard deviation. If left empty (using ""), the noise will be estimated using [Ponomarenko](https://www.ipol.im/pub/art/2013/45/).
4. output estimated kernel file. Should be a tif file in order to preserve floating point values.
5. output denoised image file.
6. output deblurring file.

Compilation
===========

This project requires some system libraries, such as `libtiff`, `libjpeg`, `libpng`, `libfftw3` and their respective headers. A recent C++ compiler, `make` and `cmake` are also required.
The programs use [iio](https://github.com/mnhrdt/iio) to read and write images.

To compile the C++ parts, use the following commands:
```bash
$ make -C ponomarenko-noise-estimation/
$ make -C kernel-estimation/
$ mkdir nonblind-deconvolution/build
$ cd nonblind-deconvolution/build
$ cmake ..
$ make
```

For FFDNet, you may use use conda and install the dependencies as such: (maybe)
```bash
$ # optional: conda create -n <env> python=3
$ conda activate <env>
$ conda install pytorch=0.4.1 -c pytorch
$ pip install -r ffdnet-denoiser/requirements.txt
```

References
==========

If you find this work useful, please consider citing it as:
> Anger, Jérémy, Mauricio Delbracio, and Gabriele Facciolo. "Efficient Blind Deblurring under High Noise Levels.", Published at International Symposium on Image and Signal Processing and Analysis (ISPA 2019).

Also available on [arxiv](https://arxiv.org/abs/1904.09154).

This method is based on the L0 kernel estimation family of methods.
For an overview and simple implementation, please refer to [Blind Image Deblurring using the l0 Gradient Prior](https://www.ipol.im/pub/art/2019/243/).

