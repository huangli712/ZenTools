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

# Get an uniform 𝑘-mesh
kmesh, weight = qeio_kmesh("dft")

# Determine the fermi level
fermi = qeio_fermi("dft", false)

# Get tight-binding hamiltonian H(𝑟)
rdeg, rvec, hamr = w90_read_hamr("dft")

# Build H(𝑘) along high-symmetry directions
println("Generate H(𝑘) where 𝑘 along high-symmetry directions")
hamk = w90_make_hamk(kpath, rdeg, rvec, hamr)

# Calculate the band structures
println("Diagonalize H(𝑘) where 𝑘 along high-symmetry directions")
eigs, evec = w90_diag_hamk(hamk)

# Build H(𝑘) in an uniform 𝑘-mesh
println("Generate H(𝑘) in an uniform 𝑘-mesh")
hamk = w90_make_hamk(kmesh, rdeg, rvec, hamr)

# Perform 𝑘-summation to calculate band levels
println("Try to evaluate the band levels")
nwann, _, nkpt = size(hamk)
level = zeros(C64, nwann)
for k = 1:nkpt
    for b = 1:nwann
        level[b] = level[b] + hamk[b,b,k] * weight[k]
    end
end
@. level = level / nkpt
@show level .- fermi

# Dump the band structures
println("Dump band structures into band.dat")
open("band.dat", "w") do fout
    nband, nkpt = size(eigs)
    for b = 1:nband
        for k = 1:nkpt
            println(fout, xpath[k], " ", eigs[b,k] - fermi)
        end
        println(fout)
    end
end

# Dump the 𝑘-list
open("kpath.dat", "w") do fout
    nband, nkpt = size(eigs)
    for k = 1:nkpt
        @printf(fout, "%12.6f %8.6f %8.6f %6.4f\n", kpath[k,:]..., 1.00)
    end
end
