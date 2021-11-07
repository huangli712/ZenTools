#!/usr/bin/env julia

#
# Project : Begonia
# Source  : band4vasp.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/07
#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

# Read the Kohn-Sham dataset
D = ir_read("dft")

# Check the validity of the `D` dict
key_list = [:MAP, :PG, :PW,
            :latt, :kmesh, :weight,
            :enk, :occupy, :Fchipsi, :fermi]
for k in key_list
    @assert haskey(D, k)
end

# Generate 𝑟-points in Wigner-Seitz cell
rdeg, rvec = w90_make_rcell(D[:latt])

# Build high-symmetry 𝑘-path
kstart = [0.0 0.0 0.0; # Γ
          0.5 0.0 0.0; # X
          0.5 0.5 0.0; # M
          0.0 0.0 0.0] # Γ
kend   = [0.5 0.0 0.0; # X
          0.5 0.5 0.0; # M
          0.0 0.0 0.0; # Γ
          0.5 0.5 0.5] # R
kpath, xpath = w90_make_kpath(100, kstart, kend)

# Build the hamiltonian in an uniform 𝑘-mesh
hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])

# Get nspin
_, _, nspin = size(D[:enk])

for p in eachindex(hamk)
    for s = 1:nspin
        # Calculate H(𝑟)
        HR = w90_make_hamr(D[:kmesh], rvec, hamk[p][:,:,:,s])

        # Build H(𝑘) along high-symmetry directions
        HK = w90_make_hamk(kpath, rdeg, rvec, HR)

        # Calculate and output the band structures
        eigs, evec = w90_diag_hamk(HK)
        nband, nkpt = size(eigs)

        # Dump the band structures
        open("newtest.dat", "w") do fout
            for b = 1:nband
                for k = 1:nkpt
                    println(fout, xpath[k], " ", eigs[b,k])
                end
                println(fout)
            end
        end # END OF IOSTREAM
    end # END OF S LOOP
end # END OF P LOOP
