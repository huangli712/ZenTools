    # Are the projectors correct?
    #
    # We will try to calculate some physical quantitites, which
    # will be written to external files or terminal for reference.
    #
    # These physical quantities include density matrix, overlap
    # matrix, local hamiltonian, full hamiltonian, and partial
    # density of states. Of course, it is time-comsuming to do
    # these things. So it is a good idea to turn off this feature
    # if everything is on the way.
    isinteractive() &&
    isfile(query_case()*".test") &&
    plo_monitor(D)

    # Are the projectors correct?
    #
    # We will try to calculate some physical quantitites, which
    # will be written to external files or terminal for reference.
    isinteractive() &&
    isfile(query_case()*".test") &&
    wannier_monitor(D)


"""
    plo_monitor(D::Dict{Symbol,Any})

Generate some key physical quantities by using the projectors and the
Kohn-Sham band structures. It is used for debug only.

See also: [`plo_adaptor`](@ref).
"""
function plo_monitor(D::Dict{Symbol,Any})
    if haskey(D, :MAP)
        # If D[:MAP] is ready, it means that D[:PW] is created and the
        # projectors are normalized and orthogonalized.

        # Calculate and output overlap matrix
        ovlp = calc_ovlp(D[:PW], D[:Fchipsi], D[:weight])
        view_ovlp(D[:PG], ovlp)

        # Calculate and output density matrix
        dm = calc_dm(D[:PW], D[:Fchipsi], D[:weight], D[:occupy])
        view_dm(D[:PG], dm)

        # Calculate and output effective atomic level
        level = calc_level(D[:PW], D[:Fchipsi], D[:weight], D[:enk])
        view_level(D[:PG], level)

        # Calculate and output full hamiltonian
        hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])
        view_hamk(hamk)

        # Calculate and output density of states
        if get_d("smear") == "tetra"
            mesh, dos = calc_dos(D[:PW], D[:Fchipsi], D[:itet], D[:enk])
            view_dos(mesh, dos)
        end
    else
        # If D[:MAP] is not ready, it means that the projectors have not
        # been postprocessed.

        # Calculate and output overlap matrix
        ovlp = calc_ovlp(D[:chipsi], D[:weight])
        view_ovlp(ovlp)

        # Calculate and output density matrix
        dm = calc_dm(D[:chipsi], D[:weight], D[:occupy])
        view_dm(dm)
    end
end

"""
    wannier_monitor(D::Dict{Symbol,Any})

Try to check and examine whether the obtained wannier functions are
correct and reasonable. Be careful, the calc_ovlp(), calc_dm(), and
calc_level() functions are defined in plo.jl.

See also: [`wannier_adaptor`](@ref).
"""
function wannier_monitor(D::Dict{Symbol,Any})
    # Calculate and output overlap matrix
    ovlp = calc_ovlp(D[:PW], D[:Fchipsi], D[:weight])
    view_ovlp(D[:PG], ovlp)

    # Calculate and output density matrix
    dm = calc_dm(D[:PW], D[:Fchipsi], D[:weight], D[:occupy])
    view_dm(D[:PG], dm)

    # Calculate and output effective atomic level
    level = calc_level(D[:PW], D[:Fchipsi], D[:weight], D[:enk])
    view_level(D[:PG], level)

    # Calculate and output full hamiltonian
    hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])
    view_hamk(hamk)
end
