#!/usr/bin/env julia

#
# Project : Begonia
# Source  : onlydiag.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/08
#

#=
*Remarks*:

This script is used to parse the matrix functions, which are essential
output of the dmft engine (`Dyson`), and extract the diagonal elements
only. Only for debug purpose.
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

cfg = inp_toml("../SrVO3.toml", true)
rev_dict(cfg)
ai = GetImpurity()

fmesh, Delta = read_delta(ai, "dmft.green")
_, qdim, nmesh, nspin, nsite = size(Delta)

file = "dmft.green.diag"
open(file, "w") do fout
    for m = 1:nmesh
        @printf(fout, "%6i%16.8f", m, fmesh[m])
        for q = 1:qdim
            z = Delta[q,q,m,1,1]
            @printf(fout, "%16.8f%16.8f", real(z), imag(z))
        end
        println(fout)
    end
end
