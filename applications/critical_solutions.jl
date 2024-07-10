# findCriticalSolution(N) takes an int starting from 1 to find the symmetric
#   piecewise polynomial having total length 2*(N+1) that is N-times continuously
#   differentiable and has maximum rate of stop-band attenuation 1/w^(2N+1)

function findCriticalSolution(N)

    # there are two elements per value of A, the value of k is the same
    k_vec1 = zeros(Rational{BigInt}, N+1)
    k_vec1[N] = 1
    k_vec2 = zeros(Rational{BigInt}, N+1)
    k_vec2[N+1] = 1
    elements = Vector{element}()
    for a in 1:(N+1)
        d1 = B(Rational{BigInt}, N, a)*(J(Rational{BigInt}, N, a)\k_vec1)
        el1 = u(Rational{BigInt}, N, 'e', 0, a, d1)
        d2 = B(Rational{BigInt}, N, a)*(J(Rational{BigInt}, N, a)\k_vec2)
        el2 = u(Rational{BigInt}, N, 'e', 0, a, d2)
        push!(elements, deepcopy(el1))
        push!(elements, deepcopy(el2))
    end

    # Matrix of constraints is 2*(N+1) by 2*(N+1)
    # standard i,j used for indexing: i is the row location, j is the column
    # rhs vector has length 2*(N+1)

    rhs = zeros(Rational{BigInt}, 2*(N+1))
    constraints = zeros(Rational{BigInt}, 2*(N+1), 2*(N+1))

    # first constraint is the sum of values at t=0 is 1
    rhs[1]=1
    for j in 1:(2*(N+1))
        constraints[1,j] = evaluate(elements[j], 0)
    end

    # the next N rows set constraints that the sum of values is 0 at t=1,...,N
    # rhs values are already 0
    for i in 1:N
        for j in (2*i+1):(2*(N+1))
            constraints[i+1,j] = evaluate(elements[j], i)
        end
    end

    # the (N+2)th constraint sets lim w->0 (sum Spectrums) -> 1

    rhs[N+2] = 1
    for j in 1:(2*(N+1))
        constraints[N+2, j] = evalSpectrum(elements[j], 0)
    end

    # the next (N-1) constraints ensure that the odd-degree derivatives are 0
    # note indexing starts from 1, u[1] is the 0th derivatives
    # u[2] is the 1st derivative
    # rhs values are already 0
    # if N=1, skip to last constraint

    row = N+3
    if N != 1
        for M in 2:2:(N+1)
            for j in 1:(2*(N+1))
                constraints[row, j] = elements[j].u[M]
            end
            row += 1
        end

        # higher order derivatives are held in the h vector
        first = iseven(N) ? 1 : 2
        for M in first:2:(N-2)
            for j in 1:(2*(N+1))
                constraints[row, j] = elements[j].h[M]
            end
            row += 1
        end
    end

    # the last constraint sets the value of the derivative at t=1
    rhs[row] = -1 * N//(N+1)
    for j in 3:(2*(N+1))
        constraints[row, j] = evalDer(elements[j], 1)
    end

    coeff = constraints\rhs

    result = deepcopy(elements)
    for i in 1:(2*(N+1))
        scalarMultiplyByReference(result[i], coeff[i])
    end

    return sumOfElements(result), elements, constraints, rhs

end
