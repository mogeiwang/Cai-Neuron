# !/usr/bin/env bash
# -*- coding:utf-8 -*-

echo "gcc..."
gcc main.c -o nif -g -lm -lGL -lGLU -lglut -Wall

echo "begin..."
./nif 10 4    1      1024000    4  >  output_run_sh.txt
gnuplot plot.gp
okular plot1-4.pdf &

echo "group -"
#     RK2/4 adjust  clicks    m
./nif 10 4    1      1000    4  >>  output_run_sh.txt
./nif 10 4    1      2000    4  >>  output_run_sh.txt
./nif 10 4    1      4000    4  >>  output_run_sh.txt
./nif 10 4    1      8000    4  >>  output_run_sh.txt
./nif 10 4    1      16000   4  >>  output_run_sh.txt
./nif 10 4    1      32000   4  >>  output_run_sh.txt
./nif 10 4    1      64000   4  >>  output_run_sh.txt
./nif 10 4    1      128000  4  >>  output_run_sh.txt

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
