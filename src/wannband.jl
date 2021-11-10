#!/usr/bin/env julia

#
# Project : Begonia
# Source  : wannband.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/10
#

#=
*Remarks*:

This script is used to analyze the projectors generated by using the
maximally localized wannier function scheme. At first, it will try
to reproduce the band structures along selected high-symmetry directions
via the wannier interpolation. Then it will evaluate the band levels.

This script is only for debug purpose. Perhaps it is only suitable for
the quantum espresso + wannier90 mode.
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

"""
    build_uniform_kmesh(x::SpecialPointsCard)

Try to generate an uniform 𝑘-mesh via SpecialPointsCard. If you can not
access regular 𝑘-mesh from the standout output of the DFT engine, perhaps
you can try this function. Note that the SpecialPointsCard struct has
been announced at ZenCore/qe.jl.
"""
function build_uniform_kmesh(x::SpecialPointsCard)
    # Print the header
    println("Generate an uniform 𝑘-mesh")

    nkpt = length(x.data)
    kmesh = zeros(F64, nkpt, 3)
    weight = zeros(F64, nkpt)
    #
    for k = 1:nkpt
        kmesh[k,:] .= x.data[k].coord
        weight[k] = x.data[k].weight
    end

    # Print some useful information
    println("  > Number of 𝑘-points: ", nkpt)
    println("  > Shape of Array kmesh: ", size(kmesh))
    println("  > Shape of Array weight: ", size(weight))

    # Return the desired arrays
    return kmesh, weight
end

"""
    calc_band_level(hamk::Array{C64,3}, weight::Array{F64,1})

Try to calculate band levels via 𝑘-summation.
"""
function calc_band_level(hamk::Array{C64,3}, weight::Array{F64,1})
    # Print the header
    println("Compute the band levels")

    nband, _, nkpt = size(hamk)
    level = zeros(C64, nband)
    #
    for k = 1:nkpt
        for b = 1:nband
            level[b] = level[b] + hamk[b,b,k] * weight[k]
        end
    end
    #
    level = level / sum(weight)

    # Print some useful information
    println("  > Number of 𝑘-points: ", nkpt)
    println("  > Number of wannier bands: ", nband)

    # Return the desired array
    return level
end

# Build high-symmetry 𝑘-path
# Number of 𝑘-points per direction. You can modify it.
ndiv = 100
#
# Please modify the following 𝑘-points to define high-symmetry directions
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
kmesh, weight = build_uniform_kmesh(SpecialPointsCard(12))
#
# Alternatively, you can use the irio_kmesh() function. But you have to
# make sure the obtained 𝑘-mesh is uniform.
#kmesh, weight = irio_kmesh("dft")

# Determine the fermi level
fermi = irio_fermi("dft")

# Get tight-binding hamiltonian H(𝑟)
rdeg, rvec, hamr = w90_read_hamr("dft")

# Build H(𝑘) along high-symmetry directions
hamk = w90_make_hamk(kpath, rdeg, rvec, hamr)

# Calculate the band structures
eigs, evec = w90_diag_hamk(hamk)

# Build H(𝑘) in an uniform 𝑘-mesh
hamk = w90_make_hamk(kmesh, rdeg, rvec, hamr)

# Perform 𝑘-summation to calculate band levels
calc_band_level(hamk, weight)

# Dump the band structures
println("Dump band structures into band.dat")
open("band.dat", "w") do fout
    nband, nkpt = size(eigs)
    for b = 1:nband
        for k = 1:nkpt
            @printf(fout, "%12.6f %12.6f\n", xpath[k], eigs[b,k] - fermi)
        end
        println(fout)
    end
end

# Dump the 𝑘-list
println("Dump 𝑘-path into kpath.dat")
open("kpath.dat", "w") do fout
    nband, nkpt = size(eigs)
    for k = 1:nkpt
        @printf(fout, "%12.6f %8.6f %8.6f %6.4f\n", kpath[k,:]..., 1.00)
    end
end

# Dump the band levels
println("Dump band levels into level.dat")
open("level.dat", "w") do fout
    for i in eachindex(level)
        @printf(fout, "%4i %12.6f %12.6f\n", i, real(level[i]), real(level[i]) - fermi)
    end
end
