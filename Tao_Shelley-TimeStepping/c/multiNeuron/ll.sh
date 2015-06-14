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
# · Description: run the multi neuron program.
# ·
# ==============================================

echo "gcc..."
gcc main.c -o nif -g -lm -lGL -lGLU -lglut -Wall

echo "begin..."
./nif 10 4    1      1024000    5  >  output_run_sh.txt
gnuplot plot.gp
okular plot1-4.pdf &

echo "group 1"
#     RK2/4 adjust  clicks    m
./nif 10 2    0      1000    5  >>  output_run_sh.txt
./nif 10 2    0      2000    5  >>  output_run_sh.txt
./nif 10 2    0      4000    5  >>  output_run_sh.txt
./nif 10 2    0      8000    5  >>  output_run_sh.txt
./nif 10 2    0      16000   5  >>  output_run_sh.txt
./nif 10 2    0      32000   5  >>  output_run_sh.txt
./nif 10 2    0      64000   5  >>  output_run_sh.txt
./nif 10 2    0      128000  5  >>  output_run_sh.txt

echo "group 2"
#     RK2/4 adjust  clicks    m
./nif 10 2    1      1000    5  >>  output_run_sh.txt
./nif 10 2    1      2000    5  >>  output_run_sh.txt
./nif 10 2    1      4000    5  >>  output_run_sh.txt
./nif 10 2    1      8000    5  >>  output_run_sh.txt
./nif 10 2    1      16000   5  >>  output_run_sh.txt
./nif 10 2    1      32000   5  >>  output_run_sh.txt
./nif 10 2    1      64000   5  >>  output_run_sh.txt
./nif 10 2    1      128000  5  >>  output_run_sh.txt

echo "group 3"
#     RK2/4 adjust  clicks    m
./nif 10 4    0      1000    5  >>  output_run_sh.txt
./nif 10 4    0      2000    5  >>  output_run_sh.txt
./nif 10 4    0      4000    5  >>  output_run_sh.txt
./nif 10 4    0      8000    5  >>  output_run_sh.txt
./nif 10 4    0      16000   5  >>  output_run_sh.txt
./nif 10 4    0      32000   5  >>  output_run_sh.txt
./nif 10 4    0      64000   5  >>  output_run_sh.txt
./nif 10 4    0      128000  5  >>  output_run_sh.txt

echo "group 4"
#     RK2/4 adjust  clicks    m
./nif 10 4    1      1000    5  >>  output_run_sh.txt
./nif 10 4    1      2000    5  >>  output_run_sh.txt
./nif 10 4    1      4000    5  >>  output_run_sh.txt
./nif 10 4    1      8000    5  >>  output_run_sh.txt
./nif 10 4    1      16000   5  >>  output_run_sh.txt
./nif 10 4    1      32000   5  >>  output_run_sh.txt
./nif 10 4    1      64000   5  >>  output_run_sh.txt
./nif 10 4    1      128000  5  >>  output_run_sh.txt

echo
echo "Neuron #0...."
echo
cat output_run_sh.txt| grep "Neuron #0"
echo
echo "Neuron #1...."
echo
cat output_run_sh.txt| grep "Neuron #1"
echo
echo "Neuron #2...."
echo
cat output_run_sh.txt| grep "Neuron #2"
echo
echo "Neuron #3...."
echo
cat output_run_sh.txt| grep "Neuron #3"
