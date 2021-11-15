#!/usr/bin/env julia

#
# Project : Begonia
# Source  : ircheck.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/15
#

#=
*Remarks*:

Are the projectors correct? It is an important question. Here we will
try to calculate some physical quantitites, which will be written to
external files or terminal for reference.

These physical quantities include density matrix, overlap matrix, local
hamiltonian, full hamiltonian, and partial density of states. Of course,
it is time-comsuming to do these things.

This script is only for debug purpose.
=#

"""
    ir_check(D::Dict{Symbol,Any})

Generate some key physical quantities by using the projectors and the
Kohn-Sham band structures. It is used for debug only.
"""
function ir_check(D::Dict{Symbol,Any})
    # Calculate and output overlap matrix
    ovlp = calc_ovlp(D[:chipsi], D[:weight])
    view_ovlp(ovlp)
    #
    ovlp = calc_ovlp(D[:PW], D[:Fchipsi], D[:weight])
    view_ovlp(D[:PG], ovlp)

    # Calculate and output density matrix
    dm = calc_dm(D[:chipsi], D[:weight], D[:occupy])
    view_dm(dm)
    #
    dm = calc_dm(D[:PW], D[:Fchipsi], D[:weight], D[:occupy])
    view_dm(D[:PG], dm)

    # Calculate and output Kohn-Sham band level
    level = calc_level(D[:chipsi], D[:weight], D[:enk])
    view_level(level)
    #
    level = calc_level(D[:PW], D[:Fchipsi], D[:weight], D[:enk])
    view_level(D[:PG], level)

    # Calculate and output hamiltonian matrix in local basis
    hamk = calc_hamk(D[:chipsi], D[:enk])
    view_hamk(hamk)
    #
    hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])
    view_hamk(hamk)

    # Calculate and output density of states
    if get_d("smear") == "tetra"
        mesh, dos = calc_dos(D[:PW], D[:Fchipsi], D[:itet], D[:enk])
        view_dos(mesh, dos)
    end
end
