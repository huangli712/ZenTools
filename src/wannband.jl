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
Perhaps this script is only suitable for qe + wannier90 mode
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

# Build high-symmetry 𝑘-path
println("Generate the high-symmetry 𝑘-path in the Brillouin zone")
#
# Number of 𝑘-points per direction
ndiv = 100
#
# Please modify the following 𝑘-points to define high-symmetry 𝑘-paths
kstart = [0.0 0.0 0.0; # Γ
          0.5 0.0 0.0; # X
          0.5 0.5 0.0; # M
          0.0 0.0 0.0] # Γ
kend   = [0.5 0.0 0.0; # X
          0.5 0.5 0.0; # M
          0.0 0.0 0.0; # Γ
          0.5 0.5 0.5] # R
#
# Generate 𝑘-list
kpath, xpath = w90_make_kpath(ndiv, kstart, kend)

fermi = qeio_fermi("dft", false)
rdeg, rvec, hamr = w90_read_hamr("dft")
sorry()

    hamk = w90_make_hamk(kpath, rdeg, rvec, hamr)

    eigs, evec = w90_diag_hamk(hamk)
    nband, nkpt = size(eigs)
    open("test.dat", "w") do fout
        for b = 1:nband
            for k = 1:nkpt
                println(fout, xpath[k], " ", eigs[b,k] - fermi)
            end
            println(fout)
        end
    end

#=
    level = zeros(C64, nband, nband)
    for k = 1:nkpt
        @. level = level + hamk[:,:,k]
    end
    @. level = level / nkpt
    for b = 1:nband
        @show b, level[b,b] - fermi
    end
=#

    enk = qeio_band("dft")
    nband, nkpt, _ = size(enk)
    @assert length(xpath) == nkpt
    open("bands.dat", "w") do fout
        for b = 1:nband
            for k = 1:nkpt
                println(fout, xpath[k], " ", enk[b,k,1] - fermi)
            end
            println(fout)
        end
    end

#=
    open("kpath.dat", "w") do fout
        for k = 1:nkpt
            @printf(fout, "%12.6f %8.6f %8.6f %6.4f\n", kpath[k,:]..., 1.00)
        end
    end
=#
