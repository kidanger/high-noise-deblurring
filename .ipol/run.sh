#!/usr/bin/env bash

set -eu

img=$1
k=$2
den=$3
deb=$4
add_noise=$5
sigma=$6

cp $img noisy.png
if [ "$add_noise" = "True" ]; then
	echo Adding noise...
	awgn $sigma noisy.png noisy.png
fi

time bash $bin/run.sh noisy.png "" "" kernel.tif $den $deb

upsa 8 0 kernel.tif - | qauto - $k

