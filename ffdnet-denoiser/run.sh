#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

input=$1
sigma=$2
output=$3

python3 -W ignore $DIR/test_ffdnet_ipol.py --input $input --noise_sigma $sigma --output $output --no_gpu

