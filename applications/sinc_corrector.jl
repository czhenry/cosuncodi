using Integrals

function f1(n, t)
    big(pi)^float(n-1)*( (n%4 == 1) ? cos(big(pi)*t) : ((n%4 == 2) ? -1*sin(big(pi)*t) : ((n%4 == 3) ? -1*cos(big(pi)*t) : sin(big(pi)*t)    ) ) )
end

function sincNDer(n, t)
    if t == 0
        return 1.0/big(n+1)*f1(n+1, 0)
    else
        if n == 0
            return sinc(big(t))
        else
            return (1.0/big(t))*(f1(n,t) - n*sincNDer(n-1, t))
        end
    end
end

function findSincCorrector(N)

    v_vec = zeros(BigFloat, 2*N)

    for i in 0:(2*N-1)
        v_vec[i+1] = -1*sincNDer(i, N+1)
    end

    velement = v(BigFloat, 2*N-1, 'e', 0, N+1, v_vec)


        # there are two elements per value of A, the value of k is the same
    k_vec1 = zeros(BigFloat, N+1)
    k_vec1[N] = 1
    k_vec2 = zeros(BigFloat, N+1)
    k_vec2[N+1] = 1
    elements = Vector{element}()
    for a in 1:(N+1)
        d1 = B(BigFloat, N, a)*(J(BigFloat, N, a)\k_vec1)
        el1 = u(BigFloat, N, 'e', 0, a, d1)
        d2 = B(BigFloat, N, a)*(J(BigFloat, N, a)\k_vec2)
        el2 = u(BigFloat, N, 'e', 0, a, d2)
        push!(elements, deepcopy(el1))
        push!(elements, deepcopy(el2))
    end

    # Matrix of constraints is 2*(N+1) by 2*(N+1)
    # standard i,j used for indexing: i is the row location, j is the column
    # rhs vector has length 2*(N+1)

    rhs = zeros(BigFloat, 2*(N+1))
    constraints = zeros(BigFloat, 2*(N+1), 2*(N+1))

    # first constraint is the sum of values at t=0 is 0
    rhs[1]=0
    for j in 1:(2*(N+1))
        constraints[1,j] = evaluate(elements[j], 0)
    end

    # the next N rows set constraints that the sum of values is 0 at t=1,...,N
    # rhs values are -1*evaluate(velement, t)
    for i in 1:N
        for j in (2*i+1):(2*(N+1))
            constraints[i+1,j] = evaluate(elements[j], i)
        end
        rhs[i+1] = -1*evaluate(velement, i)
    end

    # the (N+2)th constraint sets lim w->0 (sum Spectrums) -> 1
    # we need a numerical integral of sinc(t) on [-2(N+1), 2(N+1)]
    g(t, p) = sinc(big(t))
    domain = (-1*(N+1), (N+1))
    problem = IntegralProblem(g, domain)
    sinc_integral = solve(problem, HCubatureJL(); reltol = 1e-15, abstol = 1e-15)

    rhs[N+2] = 1 - sinc_integral.u - evalSpectrum(velement, 0)
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
    #rhs[row] = float(-1 * big(N//(N+1))) - sincNDer(1, 1)
    rhs[row] = float(-1 * big(2*N//(2*N+1))) - sincNDer(1, 1)
    for j in 3:(2*(N+1))
        constraints[row, j] = evalDer(elements[j], 1)
    end


    coeff = constraints\rhs

    result = deepcopy(elements)
    for i in 1:(2*(N+1))
        scalarMultiplyByReference(result[i], coeff[i])
    end

    push!(result, velement)

    return sumOfElements(result), elements, constraints, rhs

end
