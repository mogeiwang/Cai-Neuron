# !/usr/bin/env bash
# -*- coding:utf-8 -*-

# ==============================================
# ·
# · Author: Mogei Wang
# ·
# · MogeiWang@GMail.com
# ·
# · Filename: run.sh
# ·
# · COPYRIGHT 2015
# ·
# · Description: run the single neuron program.
# ·
# ==============================================

echo "begin..."
echo

./singleNeuron > output_run_sh.txt

#     RK2/4 adjust  clicks
./singleNeuron 2    0      1000  >> output_run_sh.txt
./singleNeuron 2    0      2000  >> output_run_sh.txt
./singleNeuron 2    0      4000  >> output_run_sh.txt
./singleNeuron 2    0      8000  >> output_run_sh.txt
./singleNeuron 2    0      16000  >> output_run_sh.txt
./singleNeuron 2    0      32000  >> output_run_sh.txt
./singleNeuron 2    0      64000  >> output_run_sh.txt
./singleNeuron 2    0      128000  >> output_run_sh.txt

#     RK2/4 adjust  clicks
./singleNeuron 2    1      1000  >> output_run_sh.txt
./singleNeuron 2    1      2000  >> output_run_sh.txt
./singleNeuron 2    1      4000  >> output_run_sh.txt
./singleNeuron 2    1      8000  >> output_run_sh.txt
./singleNeuron 2    1      16000  >> output_run_sh.txt
./singleNeuron 2    1      32000  >> output_run_sh.txt
./singleNeuron 2    1      64000  >> output_run_sh.txt
./singleNeuron 2    1      128000  >> output_run_sh.txt

#     RK2/4 adjust  clicks
./singleNeuron 4    0      1000  >> output_run_sh.txt
./singleNeuron 4    0      2000  >> output_run_sh.txt
./singleNeuron 4    0      4000  >> output_run_sh.txt
./singleNeuron 4    0      8000  >> output_run_sh.txt
./singleNeuron 4    0      16000  >> output_run_sh.txt
./singleNeuron 4    0      32000  >> output_run_sh.txt
./singleNeuron 4    0      64000  >> output_run_sh.txt
./singleNeuron 4    0      128000  >> output_run_sh.txt

#     RK2/4 adjust  clicks
./singleNeuron 4    1      1000  >> output_run_sh.txt
./singleNeuron 4    1      2000  >> output_run_sh.txt
./singleNeuron 4    1      4000  >> output_run_sh.txt
./singleNeuron 4    1      8000  >> output_run_sh.txt
./singleNeuron 4    1      16000  >> output_run_sh.txt
./singleNeuron 4    1      32000  >> output_run_sh.txt
./singleNeuron 4    1      64000  >> output_run_sh.txt
./singleNeuron 4    1      128000 >> output_run_sh.txt

echo
echo "begin..."
