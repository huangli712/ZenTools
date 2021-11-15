#!/usr/bin/env julia

#
# Project : Begonia
# Source  : ploband.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/15
#

#=
*Remarks*:
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

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

# Read the Kohn-Sham dataset
D = ir_read("dft")

# Check the validity of the `D` dict
key_list = [:MAP, :PG, :PW,
            :latt, :kmesh, :weight,
            :enk, :occupy, :Fchipsi, :fermi, :chipsi]
for k in key_list
    @assert haskey(D, k)
end

# Generate 𝑟-points in Wigner-Seitz cell
rdeg, rvec = w90_make_rcell(D[:latt])

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

# Build the hamiltonian in an uniform 𝑘-mesh
hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])

# Get nspin
_, _, nspin = size(D[:enk])

# Go through the groups and spins
for p in eachindex(hamk)
    for s = 1:nspin
        println(repeat("=", 20))
        println("Group: [$p] Spin: [$s]")
        println(repeat("=", 20))

        # Calculate H(𝑟)
        HR = w90_make_hamr(D[:kmesh], rvec, hamk[p][:,:,:,s])

        # Build H(𝑘) along high-symmetry directions
        HK = w90_make_hamk(kpath, rdeg, rvec, HR)

        # Calculate the band structures
        eigs, evec = w90_diag_hamk(HK)
        nband, nkpt = size(eigs)

        # Dump the band structures
        println("  > Dump band structures into band.plo.p$p.s$s")
        open("band.plo.p$p.s$s", "w") do fout
            for b = 1:nband
                for k = 1:nkpt
                    println(fout, xpath[k], " ", eigs[b,k])
                end
                println(fout)
            end
        end # END OF IOSTREAM
    end # END OF S LOOP
end # END OF P LOOP

# Build the hamiltonian in an uniform 𝑘-mesh
hamk = calc_hamk(D[:chipsi], D[:enk])

# Go through the spins
for s = 1:nspin
    println(repeat("=", 20))
    println("Spin: [$s]")
    println(repeat("=", 20))

    # Perform 𝑘-summation to calculate band levels
    level = calc_band_level(hamk[:,:,:,s], D[:weight])

    # Calculate H(𝑟)
    HR = w90_make_hamr(D[:kmesh], rvec, hamk[:,:,:,s])

    # Build H(𝑘) along high-symmetry directions
    HK = w90_make_hamk(kpath, rdeg, rvec, HR)

    # Calculate the band structures
    eigs, evec = w90_diag_hamk(HK)
    nband, nkpt = size(eigs)

    # Dump the band structures
    println("  > Dump band structures into band.plo.s$s")
    open("band.plo.s$s", "w") do fout
        for b = 1:nband
            for k = 1:nkpt
                println(fout, xpath[k], " ", eigs[b,k])
            end
            println(fout)
        end
    end # END OF IOSTREAM

    # Dump the band levels
    println("Dump band levels into level.plo.s$s")
    open("level.plo.s$s", "w") do fout
        for i in eachindex(level)
            @printf(fout, "%4i %12.6f\n", i, real(level[i]))
        end
    end # END OF IOSTREAM
end # END OF S LOOP

# Dump the 𝑘-list
println("Dump 𝑘-path into kpath.plo")
open("kpath.plo", "w") do fout
    nkpt, _ = size(kpath)
    for k = 1:nkpt
        @printf(fout, "%12.6f %8.6f %8.6f %6.4f\n", kpath[k,:]..., 1.00)
    end
end # END OF IOSTREAM
