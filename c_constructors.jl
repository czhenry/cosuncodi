export u, v, w, nullElement, uBasis, vBasis

#add function to be used in natBasis constructor and evalDer on left/right



# functions u, v, and w are constructors of elements
# these functions match with the notation
# constructing an element with u, v, or w fills in c, c_hat, h, k, n, and spec

function u(numType, N::Int64, sym::Char, shift, a, d)
    result=element(deepcopy(numType), deepcopy(N), deepcopy(sym), deepcopy(shift), deepcopy(a), deepcopy(d), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, 2*N+2), deepcopy(nullSpectrum))
    if (a != 0)
      result.c=B(numType, N, a)\result.u
      result.h=G(numType, N, a)*result.c
      result.k=J(numType, N, a)*result.c
      result.n=M(numType, N, a)*result.c
      constructSpectrum(result)
    end
    result
end

function v(numType, N::Int64, sym::Char, shift, a, d)
    result=element(deepcopy(numType), deepcopy(N), deepcopy(sym), deepcopy(shift), deepcopy(a), zeros(numType, N+1), deepcopy(d), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, 2*N+2), deepcopy(nullSpectrum))
    if (a != 0)
      result.c_hat=B_hat(numType, N, a)\result.v
      result.h=G_hat(numType, N, a)*result.c_hat
      result.k=J_hat(numType, N, a)*result.c_hat
      result.n=M_hat(numType, N, a)*result.c_hat
      constructSpectrum(result)
    end
    result
end

function w(numType, N::Int64, sym::Char, shift, a, d1, d2)
    result=element(deepcopy(numType), deepcopy(N), deepcopy(sym), deepcopy(shift), deepcopy(a), deepcopy(d1), deepcopy(d2), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, N+1), zeros(numType, 2*N+2), deepcopy(nullSpectrum))
    if (a != 0)
      result.c=B(numType, N, a)\result.u
      result.h=G(numType, N, a)*result.c
      result.k=J(numType, N, a)*result.c
      result.c_hat=B_hat(numType, N, a)\result.v
      result.h+=G_hat(numType, N, a)*result.c_hat
      result.k+=J_hat(numType, N, a)*result.c_hat
      result.n=M(numType, N, a)*result.c
      result.n+=M_hat(numType, N, a)*result.c_hat
      constructSpectrum(result)
    end
    result
end

# add constructor for elements from vector n as n0, n1, ..., n(2*N+1)
# function nB(numType, N::Int64, sym::Char, shift, a, n)

nullElement=u(Int64, 0, 'n', 0, 0, [])


# 2 basis creators to create elements
function uBasis(numType, N::Int64, sym::Char, shift, a)
    result=collection([],0,0,0,deepcopy(nullSpectrum))
    for i in 0:N
        d = zeros(Int64, N+1)
        d[i+1] = 1
        push!(result.elements, u(numType, N, sym, shift, a, d))
    end
    result
end

function vBasis(numType, N::Int64, sym::Char, shift, a)
    result=collection([],0,0,0,deepcopy(nullSpectrum))
    for i in 0:N
        d = zeros(Int64, N+1)
        d[i+1] = 1
        push!(result.elements, v(numType, N, sym, shift, a, d))
    end
    result
end
