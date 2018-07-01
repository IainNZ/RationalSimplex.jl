module RationalSimplex

import LinearAlgebra: I

export simplex

"""
    simplex(c, obj, A, b, senses)

Solve the linear program

     {min/max}  dot(c, x)
    subject to  A x {>=|==|<=} b
                x >= 0,

where `x` is in the set of Rationals (so algorithm is exact).
`obj` should be either `:Min` or `:Max`.
`senses` should be a `Vector` of `Char`s: `'<','=','>'`.
"""
function simplex(c::Vector{T}, obj::Symbol, A::Matrix{T}, b::Vector{T},
                 senses::Vector{Char}) where {T <: Rational}
    # Any constraint that isn't an equality needs an auxiliary variable.
    num_auxiliary = count(sense -> sense != '=', senses)
    num_constraints, num_variables = size(A)
    num_variables_aux = num_variables + num_auxiliary
    # Extend constraint matrix and objective function, normalize objective.
    A_aux = zeros(T, num_constraints, num_variables_aux)
    A_aux[:, 1:num_variables] = A
    c_aux = zeros(T, num_variables_aux)
    c_aux[1:num_variables] = (obj == :Max) ? -c : c
    b_aux = copy(b)  # In case we modify signs.
    # Add the auxiliaries to the constraint matrix, and make sure b ≥ 0.
    offset = 1
    for (row, sense) in enumerate(senses)
        if sense == '<'
            A_aux[row, num_variables + offset] = one(T)
            offset += 1
        elseif sense == '>'
            A_aux[row, num_variables + offset] = -one(T)
            offset += 1
        end
        if b[row] < zero(T)
            A_aux[row, :] *= -one(T)
            b_aux[row] *= -one(T)
        end
    end
    # Solve problem with auxiliaries, but remove auxiliaries before returning.
    sense, x = simplex(c_aux, A_aux, b_aux)    
    return sense, x[1:num_variables]
end


"""
    simplex(c, A, b)

Solve the linear program in standard computational form

           min  dot(c, x)
    subject to  A x == b [≥ 0]
                x >= 0,

where `x` is in the set of Rationals (so algorithm is exact).

The algorithm is the "two-phase primal revised simplex method".
In the first phase auxiliaries are created to get an initial basis.
We then eliminate the auxiliaries from the basis until it consists only of
actual variables. This is the "textbook" algorithm, and shouldn't be used for
anything that matters. It doesn't exploit sparsity at all. You could use it
with floating points but it won't work for anything except the most simple
problem due to accumulated errors and comparisons with zero.
"""
function simplex(c::Vector{T}, A::Matrix{T}, b::Vector{T}) where {T<:Rational}
    @assert all(b .> zero(T))

    # Set up data structures.
    num_constraints, num_variables = size(A)
    is_basic = zeros(Bool, num_variables + num_constraints)
    basic = zeros(Int, num_constraints)  # Indices of current basis.
    Binv = Matrix{T}(I, num_constraints, num_constraints)  # Basis inverse.
    cB = ones(T, num_constraints)  # Costs of basic variables.
    x = zeros(T, num_variables + num_constraints)  # Current solution.

    # Intialize phase one of RSM by setting basis = auxiliaries.
    for aux_index in 1:num_constraints
        basic[aux_index] = num_variables + aux_index
        is_basic[num_variables + aux_index] = true
        x[num_variables + aux_index] = b[aux_index]
    end
    phase_one = true

    # Begin simplex iterations.
    status = :Unknown
    while true
        # Calculate dual solution.
        π = vec(cB' * Binv)
        # Use it to calculate the reduced costs of the variables. Don't
        # calculate for auxiliaries - they can't re-enter the basis.
        rc = (phase_one ? zero(T) : c) .- vec(π' * A)
        @. rc *= !is_basic[1:num_variables]  # In basis, so don't consider.
        # Find variable to enter basis.
        min_rc, entering = findmin(rc)
        # If none have negative reduced cost, there are no variables that can
        # enter, and we are at optimality for this phase. This may or may not
        # be optimal for the actual problem.
        if min_rc >= zero(T)
            if phase_one
                phase_one = false
                # If any auxiliary still nonzero, we couldn't find a feasible
                # basis without auxiliaries.
                if any(x[num_variables+1:end] .> zero(T))
                    return :Infeasible, x[1:num_variables]
                end
                # Otherwise, start phase two with nice feasible basis.
                for i in 1:num_constraints
                    cB[i] = basic[i] > num_variables ? zero(T) : c[basic[i]]
                end
                continue
            else  # Phase two. We can't improve further, so we are optimal.
                return :Optimal, x[1:num_variables]
            end
        end
        # Calculate how the solution will change when our new variable enters
        # the basis and increases from zero.
        BinvAs = Binv * A[:, entering]
        # Perform a "ratio test" on each variable to determine which will
        # reach zero first.
        leaving = 0
        min_ratio = zero(T)
        for j in 1:num_constraints
            if BinvAs[j] > zero(T)
                ratio = x[basic[j]] / BinvAs[j]
                if ratio < min_ratio || leaving == 0
                    min_ratio = ratio
                    leaving = j
                end
            end
        end
        # If no variable will leave basis, then we have an unbounded problem.
        if leaving == 0
            return :Unbounded, x[1:num_variables]
        end
        # Update solution.
        for j in 1:num_constraints
            x[basic[j]] -= min_ratio * BinvAs[j]
        end
        x[entering] = min_ratio
        # Update basis inverse.
        # Our tableau is [ Binv b | Binv | BinvAs ] and we doing a pivot on the
        # `leaving` row of BinvAs.
        pivot_value = BinvAs[leaving]
        for basis_row in 1:num_constraints
            basis_row == leaving && continue  # All rows except `leaving` row.
            factor = BinvAs[basis_row] / pivot_value
            for basis_col in 1:num_constraints
                Binv[basis_row, basis_col] -= factor * Binv[leaving, basis_col]
            end
        end
        # Finally, the `leaving` row.
        for basis_col in 1:num_constraints
            Binv[leaving, basis_col] /= pivot_value
        end
        # Update variable status flags.
        is_basic[basic[leaving]] = false
        is_basic[entering] = true
        cB[leaving] = phase_one ? zero(T) : c[entering]
        basic[leaving] = entering
    end
end

end  # module RationalSimplex.