#  !/usr/bin/env bash
#  -*- coding:utf-8 -*-

pandoc -f markdown+lhs slides.md -o slides.html -t dzslides -s -S --toc

