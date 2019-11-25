#!/usr/bin/env bash

set -e

blurry=$1
kernelsize=${2:-35}
sigma=$3
estimatedkernel=${4:-$blurry-kernel.tif}
denoised=${5:-$blurry-denoised.tif}
deblurred=${6:-$blurry-deblurred.png}

echo "Blurry image: $blurry"
echo "Kernel size: $kernelsize"

if [ -z "$sigma" ]; then
    echo "Noise level estimation..."
    sigma=$(./ponomarenko-noise-estimation/ponomarenko -b 3 $blurry | tail -n 2 | head -n 1 | cut -d ' ' -f 7)
    echo "Estimated sigma: $sigma/255"
else
    echo "User sigma: $sigma/255"
fi

# same ratios as in the paper
lambda=$(echo "0.5 * $sigma/255" | bc -l)
gamma=$(echo "200 * $sigma/255" | bc -l)
echo "Kernel estimation..."
echo " - set lambda-min to $lambda"
echo " - set gamma to $gamma"
./kernel-estimation/estimate-kernel $kernelsize $blurry $estimatedkernel --lambda-min $lambda --gamma $gamma
echo "Estimated kernel: $estimatedkernel"

echo "Denoising..."
bash ffdnet-denoiser/run.sh $blurry $sigma $denoised
echo "Denoised image: $denoised"

lambda2=0.005
echo "Non-blind deconvolution..."
./nonblind-deconvolution/build/deblur $denoised $estimatedkernel $deblurred $lambda2
echo "Deblurring result: $deblurred"

