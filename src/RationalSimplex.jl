module RationalSimplex

export simplex

#######################################################################
# simplex
# Solve the linear program
#  {min/max}  dot(c,x)
# subject to  A x {>=|==|<=} b
#               x >= 0
# where x is in the set of Rationals (so algorithm is exact).
#######################################################################
function simplex{T<:Rational}(
    c::Array{T, 1},     # Cost vector
    obj::Symbol,        # Objective sense - :Min or :Max
    A::Array{T, 2},     # Constraint matrix
    b::Array{T, 1},     # Resource vector
    sense::Array{Char}) # Constraint senses - '<', '=', '>'
    
    m, n     = size(A)

    true_b = copy(b)

    # Count number of auxiliaries we will need
    extra = 0
    for i = 1:m
        if sense[i] != '='
            extra += 1
        end
    end
    true_A = zeros(T,m,n+extra)
    true_A[:,1:n] = A
    true_c = zeros(T,  n+extra)
    true_c[1:n] = (obj == :Max) ? -c : c

    # Add the auxiliaries
    offset = 1
    for i = 1:m
        if sense[i] == '<'
            true_A[i,n+offset] =  one(T)
            offset += 1
        elseif sense[i] == '>'
            true_A[i,n+offset] = -one(T)
            offset += 1
        end
    end

    # Make sure right-hand-side is non-negative
    for i = 1:m
        if b[i] < zero(T)
            true_A[i,:] *= -one(T)
            true_b[i]   *= -one(T)
        end
    end

    sense, x = simplex(true_c, true_A, true_b)    
    return sense, x[1:n]
end


#######################################################################
# simplex
# Solve the linear program in standard computational form
#        min  dot(c,x)
# subject to  A x == b
#               x >= 0
# where x is in the set of Rationals (so algorithm is exact).
# b must be >= 0.
# 
# The algorithm is the two-phase primal revised simplex method.
# In the first phase auxiliaries are created which we eliminate
# until we have a basis consisting solely of actual variables.
# This is pretty much the "textbook algorithm", and shouldn't
# be used for anything that matters. It doesn't exploit sparsity
# at all. You could use it with floating points but it wouldn't 
# work for anything except the most simple problem due to accumulated
# errors and the comparisons with zero.
#######################################################################
function simplex{T<:Rational}(
    c::Array{T, 1},  # Cost vector
    A::Array{T, 2},  # Constraint matrix
    b::Array{T, 1})  # Resource vector
    
    m, n     = size(A)

    for i = 1:m
        @assert b[i] > zero(T)
    end

    is_basic = zeros(Bool,n+m)
    basic    = zeros(m)     # indices of current basis
    Binv     = eye(T,m)     # inverse of basis matrix
    cB       = ones(T,m)    # costs of basic variables
    x        = zeros(T,m+n) # current solution

    # Intialize phase 1
    for j = 1:m
        basic[j] = j+n
        is_basic[j+n] = true
        x[j+n] = b[j]
    end
    #println("x")
    #println(x)
    phase_one = true

    # Begin simplex iterations
    status = :Unknown
    while true
        #println("Iter")
        # Calculate dual solution...
        pi_T = vec(cB'*Binv)
        #dump(cB)
        #dump(Binv)
        #dump(pi_T)
        
        # ... and thus the reduced costs
        entering = 0
        for j = 1:n
            is_basic[j] && continue
            rc = (phase_one ? zero(T) : c[j]) - dot(pi_T, A[:,j])
            #println("s= ", j, " ", rc)
            if rc < zero(T)
                entering = j
                #println("picked")
                break
            end
        end

        # If we couldn't find a variable with a negative reduced cost, 
        # we terminate this phase because we are at optimality for this
        # phase - not necessarily optimal for the actual problem.
        if entering == 0
            if phase_one
                phase_one = false
                #println("Phase 1 over")
                # Check objective - if 0, we are OK
                any_artificial_x_nonzero = false
                for j = n+1:n+m
                    if x[j] > zero(T)
                        any_artificial_x_nonzero = true
                        break
                    end
                end
                if any_artificial_x_nonzero
                    # It couldn't reduce objective to 0 which is equivalent
                    # to saying a feasible basis with no artificials could
                    # not be found
                    return :Infeasible, x[1:n]
                end
                # Start again in phase 2 with our nice feasible basis
                for i = 1:m
                    cB[i] = basic[i] > n ? 0.0 : c[basic[i]]
                end
                continue
            else
                return :Optimal, x[1:n]
            end
        end

        # Calculate how the solution will change when our new
        # variable enters the basis and increases from 0
        BinvAs = Binv * A[:,entering]
        #println("BinvAs")
        #println(BinvAs)

        # Perform a "ratio test" on each variable to determine
        # which will reach 0 first
        leaving = 0
        min_ratio = zero(T)
        for j = 1:m
            if BinvAs[j] > zero(T)
                ratio = x[basic[j]] / BinvAs[j]
                if ratio < min_ratio || leaving == 0
                    min_ratio = ratio
                    leaving = j
                end
            end
        end

        # If no variable will leave basis, then we have an 
        # unbounded problem.
        #println("r= ", leaving, " ", min_ratio)
        if leaving == 0
            return :Unbounded, x[1:n]
        end

        # Now we update solution...
        for i = 1:m
            x[basic[i]] -= min_ratio * BinvAs[i]
        end
        x[entering] = min_ratio
        #println("x")
        #println(x)

        # ... and the basis inverse...
        # Our tableau is [ Binv b | Binv | BinvAs ]
        # and we doing a pivot on the leaving row of BinvAs
        pivot_value = BinvAs[leaving]
        for i = 1:m  # all rows except leaving row
            i == leaving && continue
            factor = BinvAs[i] / pivot_value
            for j = 1:m
                Binv[i, j] -= factor * Binv[leaving, j]
            end
        end
            for j = 1:m
                Binv[leaving, j] /= pivot_value
            end

        # ... and variable status flags
        is_basic[basic[leaving]] = false
        is_basic[entering] = true
        cB[leaving] = phase_one ? zero(T) : c[entering]
        basic[leaving] = entering

        #println("basic")
        #println(basic)
        #readline()
    end
end

end #module