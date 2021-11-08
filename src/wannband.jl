#!/usr/bin/env julia

#
# Project : Begonia
# Source  : wannband.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/08
#

#=
*Remarks*:

This script is used to generate the DFT band structures via the
maximally localized wannier function scheme. Only for debug purpose.
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

