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

function test_plo_hamk()
    # Read H(𝑘) in uniform 𝑘-mesh
    hamk = nothing
    open("dft/hamk.chk.1", "r") do fin
        readline(fin)
        readline(fin)
        readline(fin)

        ngroup = parse(I64, line_to_array(fin)[3])
        nproj = parse(I64, line_to_array(fin)[3])
        nkpt = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        #@assert ngroup == 1
        @assert nspin == 1
        readline(fin)

        hamk = zeros(C64, nproj, nproj, nkpt, nspin)
        for s = 1:nspin
            for k = 1:nkpt
                for q = 1:nproj
                    for p = 1:nproj
                        _re, _im = parse.(F64, line_to_array(fin)[1:2])
                        hamk[p,q,k,s] = _re + _im * im
                    end
                end
            end
        end
    end

    D = ir_read("dft")
    hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])

    # Get uniform 𝑘-mesh
    #kmesh, weight = vaspio_kmesh("dft")

    # Generate 𝑟-points in Wigner-Seitz cell
    #latt = vaspio_lattice("dft")
    rdeg, rvec = w90_make_rcell(D[:latt])

    # Calculate H(𝑟)
    hamr = w90_make_hamr(D[:kmesh], rvec, hamk[1][:,:,:,1])

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

    # Build H(𝑘) along high-symmetry directions
    newhamk = w90_make_hamk(kpath, rdeg, rvec, hamr)

    # Calculate and output the band structures
    eigs, evec = w90_diag_hamk(newhamk)
    nband, nkpt = size(eigs)
    open("newtest.dat", "w") do fout
        for b = 1:nband
            for k = 1:nkpt
                println(fout, xpath[k], " ", eigs[b,k])
            end
            println(fout)
        end
    end
end