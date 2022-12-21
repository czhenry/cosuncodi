# the spectrum struct stores information about the Fourier Transform
# of an element or a collection of elements

mutable struct spectralTerm
    wDeg              # terms appear as coefficients of 1/(w)^N
    a
    shift             # named to note factor of e^(-iw*shift)
    coeff             # all coeff are complex numbers
    eiwCoeff          # factor e^(-iaw) in non-symmetric, joins together in conj
    cosCoeff          # factor cos(aw)
    sinCoeff          # factor sin(aw)
end

nullTerm=spectralTerm(0, 0, 0, 0, 0, 0, 0)

mutable struct spectrum
    terms::Vector{spectralTerm}
end    

nullSpectrum=spectrum([])

# the element struct stores information about a polynomial of degree at most 2N+2
# on (shift, shift+a) by knowing its value and N derivatives at the left and right
# end points

# N is the number of times continuously differentiable at either the left or right
# end points that we set with vectors u and v

mutable struct element
    numType
    N::Int64
    sym::Char       # for symmetry, values 'e' for even, 'o' for odd, 'n' for none
    shift           # tau in notation
    a               # length of the interval on which element is non-zero
    u     # N+1, value of the element and its 1st N derivatives, left hand side
    v     # N+1, same, but for right hand side
    h     # N+1, higher order derivatives on the left
    k     # N+1, higher order derivatives on the right
    c     # N+1, basis in U, where there is (1-t/a)^(N+1) common factor
    c_hat # N+1, basis in V, where there is (t/a)^(N+1) common factor
    n     # 2N+2, vector for natural coefficients of polynomial of t
    spec::spectrum    
end

# data structure collection is a work in progress
#  basic thing is to be able to extract the quantities we want from it

mutable struct collection
    elements::Vector{Union{collection,element}}
    pointsOfContinuity
    N
    sum
    spec::spectrum
end

nullCollection=collection([], 0, 0, 0, deepcopy(nullSpectrum))

# the sumSpectrums function will add terms of the 2nd spectrum
# into the first one, and return it
# pass by reference on spec1, pass by value on spec2

# outside: use deepcopy on 1st arg to create pass by value on both
#  and return newly allocated spectrum

function sumSpectrums(spec1::spectrum, spec2::spectrum)
    if isempty(spec2.terms)
        return spec1
    end
    if isempty(spec1.terms)
        for term in spec2.terms
            push!(spec1.terms, deepcopy(term))
        end
        return spec1
    end
    for term in spec2.terms
        found = false
        index = 0
        for i in 1:length(spec1.terms)
            if (spec1.terms[i].wDeg == term.wDeg) && (spec1.terms[i].a == term.a) && (spec1.terms[i].shift == term.shift)
                found = true
                index=i
            end
        end
        if found
            spec1.terms[index].coeff+=term.coeff
            spec1.terms[index].eiwCoeff+=term.eiwCoeff
            spec1.terms[index].sinCoeff+=term.sinCoeff
            spec1.terms[index].cosCoeff+=term.cosCoeff
        else
            push!(spec1.terms, deepcopy(term))
        end
    end
    return spec1
end

# reduceSpectrum removes terms that are 0
# pass by reference
function reduceSpectrum(this::spectrum)
    for i in length(this.terms):-1:1
        if (this.terms[i].coeff == 0) && (this.terms[i].eiwCoeff == 0) && (this.terms[i].sinCoeff == 0) && (this.terms[i].cosCoeff == 0)
            popat!(this.terms, i)
        end
    end
end

# needed reliable method to get im^N in constructSpectrum
function imN(N)
    N>=0 ? (N%4<2 ? N%2*im + (N+1)%2 : N%2*-im - (N+1)%2) : (N%4>-2 ? N%2*im - (N-1)%2 : N%2*-im + (N-1)%2)
end

# constructSpectrum will recursively create spectrum in the data structure
# on all collections/elements  
# pass by reference on a collection or element
function constructSpectrum(this::Union{collection,element})
    if typeof(this) == collection
        if isempty(this.spec.terms)
            for el in this.elements
                    constructSpectrum(el)
                    sumSpectrums(this.spec, el.spec)
            end
            reduceSpectrum(this.spec)
        end
    elseif typeof(this) == element
        if isempty(this.spec.terms)    
            if this.sym == 'n'
                for i in 0:this.N
                    if (this.u[i+1] != 0) || (this.v[i+1] != 0)
                    # wDeg a shift coeff eiwCoeff cosCoeff sinCoeff  
                        new_term = spectralTerm(i+1, this.a, this.shift, imN(-1-i)*this.u[i+1], imN(-1-i)*this.v[i+1], 0, 0)
                        push!(this.spec.terms, new_term)
                    end
                end
                for i in (this.N+1):(2*this.N+1)
                    if (this.h[i-this.N] != 0) || (this.k[i-this.N] !=0) 
                        new_term = spectralTerm(i+1, this.a, this.shift, imN(-1-i)*this.h[i-this.N], -1*imN(-1-i)*this.k[i-this.N], 0, 0)
                        push!(this.spec.terms, new_term)
                    end
                end
            else 
                # FOR SYMMETRIC TYPES, there is a choice here about sym=self-adjoint
                # for complex valued polynomials
                # current: 'e' = U(f)+conj(U(f)) which is self-adjoint
                # 'e'=>u(t)+conj(u(-t)) , 'o'=>u(t)-conj(u(-t))
                for j in 0:this.N
                    i = j+1
                    if ((this.sym == 'e') && (j%2 == 0)) || ((this.sym == 'o') && (j%2 == 1))
                        if (imag(this.u[i]) != 0) || (real(this.v[i]) != 0) || (imag(this.v[i]) != 0)
                            new_term =  spectralTerm(i, this.a, this.shift, 2im*imN(-i)*imag(this.u[i]), 0, -2im*imN(-i)*imag(this.v[i]), 2im*imN(-i)*real(this.v[i]) )
                            push!(this.spec.terms, new_term)
                        end
                    elseif (real(this.u[i]) != 0) || (real(this.v[i]) != 0) || (imag(this.v[i]) != 0)
                            new_term =  spectralTerm(i, this.a, this.shift, 2*imN(-i)*real(this.u[i]), 0, -2*imN(-i)*real(this.v[i]), -2*imN(-i)*imag(this.v[i]))
                            push!(this.spec.terms, new_term)
                    end
                end
                for j in (this.N+1):(2*this.N+1)
                    i = j - this.N
                    if ((this.sym == 'e') && (j%2 == 0)) || ((this.sym == 'o') && (j%2 == 1))
                        if (imag(this.h[i]) != 0) || (real(this.k[i]) != 0) || (imag(this.h[i]) != 0)
                            new_term =  spectralTerm(j+1, this.a, this.shift, 2im*imN(-1-j)*imag(this.h[i]), 0, -2im*imN(-1-j)*imag(this.k[i]), 2im*imN(-1-j)*real(this.k[i]))
                            push!(this.spec.terms, new_term)
                        end
                    elseif (real(this.h[i]) != 0) || (real(this.k[i]) != 0) || (imag(this.k[i]) != 0)
                            new_term =  spectralTerm(j+1, this.a, this.shift, 2*imN(-1-j)*real(this.h[i]), 0, -2*imN(-1-j)*real(this.k[i]), -2*imN(-1-j)*imag(this.k[i]))
                            push!(this.spec.terms, new_term)
                    end
                end
        
            end    
        end
        
    end
end    

# Consider breaking out these linear algebra functions in another file
# Consider also adding dictionaries to store re-used matrices

function permutation(N, M)
    factorial(N)//factorial(N-M)
end

function B(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in 0:N; for j in 0:N
      if (i<j)
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-j)*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))

        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function B_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in 0:N; for j in 0:N
      if (i<j)
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^j*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))        
        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function G(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-j)*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^(i-j)*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function G_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-N-1)*binomial(big(i),big(i- N-1))*permutation(big(j), big(i-N-1))*factorial(big(N+1))        
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-N-1)*binomial(i,i- N-1)*permutation(j, i-N-1)*factorial(N+1)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^(i-N-1)*binomial(i, i-N-1)*permutation(j, i-N-1)*factorial(N+1)
      end
    end
    end
    result
end

function J(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(N+1)*binomial(big(i), big(N+1))*permutation(big(j), big(i-N-1))*factorial(big(N+1))        
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^(N+1)*binomial(i, N+1)*permutation(j, i-N-1)*factorial(N+1)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^(N+1)*binomial(i, N+1)*permutation(j, i-N-1)*factorial(N+1)
      end
    end
    end
    result
end

function J_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, N+1, N+1)
    for i in (N+1):(2*N+1); for j in 0:N
      if (i>(j+N+1))
        result[i-N,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i-N,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^j*binomial(big(i), big(j))*permutation(big(N+1), big(i-j))*factorial(big(j))        
        else
            result[i-N,j+1]=(1/rationalize(float(a))^i)*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
        end
      else
        result[i-N,j+1]=1/a^i*(-1)^j*binomial(i, j)*permutation(N+1, i-j)*factorial(j)
      end
    end
    end
    result
end

function M(numType, N::Int64, a)
    result=Matrix{numType}(undef, 2*N+2, N+1)
    for i in 0:(2*N+1); for j in 0:N
      if (i>(j+N+1))||(j>i)
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-j)*binomial(big(N+1), big(i-j))        
        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-j)*binomial(N+1, i-j)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^(i-j)*binomial(N+1, i-j)
      end
    end
    end
    result
end

function M_hat(numType, N::Int64, a)
    result=Matrix{numType}(undef, 2*N+2, N+1)
    for i in 0:(2*N+1); for j in 0:N
      if (i>(j+N+1))||(i<(N+1))
        result[i+1,j+1]=0
      elseif (numType <: Rational)
        if numType == Rational{BigInt}
            result[i+1,j+1]=(1/big(rationalize(float(a)))^i)*(-1)^(i-N-1)*binomial(big(j), big(i-N-1))        
        else
            result[i+1,j+1]=(1/rationalize(float(a))^i)*(-1)^(i-N-1)*binomial(j, i-N-1)
        end
      else
        result[i+1,j+1]=1/a^i*(-1)^(i-N-1)*binomial(j, i-N-1)
      end
    end
    end
    result
end

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

function complexToString(z)
    result=""
    if typeof(real(z)) <: Rational
        if real(z).den == 1
            strez = string(real(z).num)
        else
            strez = "\\frac{"*string(real(z).num)*"}{"*string(real(z).den)*"}"
        end
        if imag(z).den == 1
            stabimz = string(abs(imag(z).num))*"i"
        else
            stabimz = "\\frac{"*string(abs(imag(z).num))*"}{"*string(imag(z).den)*"}i"
        end
    else
        strez = string(real(z))
        stabimz = string(abs(imag(z)))*"i"
    end
    if imag(z) < 0
        if (real(z) != 0) && (imag(z) != 0)
            result *= "("
        end
        if real(z) != 0
            result *= strez
        end
        result *= "-" 
        if imag(z) != -1
            result *= stabimz
        else
            result *= "i"
        end
        if (real(z) != 0) && (imag(z) != 0)
            result *= ")"
        end
    elseif (real(z) != 0) && (imag(z) != 0)
        result *= "("
        result *= strez 
        result *= "+" 
        if imag(z) != 1
            result *= stabimz
        else
            result *= "i"
        end
        result *= ")"
    elseif (real(z) == 0) && (imag(z) == 0) 
        result *= string(0)
    else
        if real(z) != 0
            result *= strez
        end
        if imag(z) != 0
            if imag(z) != 1
                result *= stabimz
            else
                result *= "i"
            end
        end
    end
    return result
end

function strSpectrum(this::Union{collection,element})
    result = ""
    if isempty(this.spec.terms)
        constructSpectrum(this)
    end
    first_term = true
    for term in this.spec.terms
        if !first_term
            result *= "+"
        end
        first = true
        result *= "\\frac{"
        if term.coeff != 0
            result *= complexToString(term.coeff)
            first = false
        end
        if term.eiwCoeff != 0
            if !first
                result *= "+"
            end
            result *= complexToString(term.eiwCoeff)*"e^{-i"*string(term.a)*"\\omega}"
            first = false
        end
        if term.cosCoeff != 0
            if !first
                result *= "+"
            end
            result *= complexToString(term.cosCoeff)*"cos("*string(term.a)*"\\omega)"
            first = false
        end
        if term.sinCoeff != 0
            if !first
                result *= "+"
            end
            result *= complexToString(term.sinCoeff)*"sin("*string(term.a)*"\\omega)"
            first = false
        end
        result *= "}{\\omega^{"*string(term.wDeg)*"}}"
        if term.shift != 0
            result *= "*e^{-i"*string(term.shift)*"\\omega}"
        end
        first_term = false
    end
    return result
end

# Write Print functions to produce:
# 1. code that calculates these functions simply and writes to a table (Octave, C)
# 3. on screen in julia

function vectorToString(vec)
    result = "["
    for i in 1:(length(vec)-1)
        result *= complexToString(vec[i]) * ", "
    end
    result *= complexToString(vec[length(vec)]) * "]"
    return result
end

function strElement(this::Union{collection,element})
    result = ""
    if typeof(this) == element
        result *= "\\prescript^{"
        if this.shift != 0
            result *= string(this.shift)
        end
        result *= "}_{"*string(this.a)*"}{} "
        thistype = 0
        for val in this.u
            if val != 0
                thistype |= 1
            end
        end
        for val in this.v
            if val != 0
                thistype |= 2
            end
        end
        if thistype == 1
            result *= "u"
        elseif thistype == 2
            result *= "v"
        elseif thistype == 3
            result *= "w"
        end
        result *= "_{"
        if thistype == 1
            result *= vectorToString(this.u)
        elseif thistype == 2
            result *= vectorToString(this.v)
        elseif thistype == 3
            result *= vectorToString(this.u) * ", " * vectorToString(this.v)
        end
        result *= "}"
        if this.sym != 'n'
            result *= "^{"*this.sym*"}"
        end
    else
        first = true
        for el in this.elements
            if !first
                result *= "+"
            end
            result *= strElement(el)
            first = false
        end
    end
    return result
end

# add argument to evaluate from left or from right (when at end points)

function evaluate(this::Union{collection,element}, t)
    if typeof(this) == element
        if this.sym == 'n'
            if (this.shift <= t <= this.shift+this.a)
                return evalpoly(t-this.shift, this.n)
            else
                return 0
            end
        else
            # We must be symmetric
            if (this.shift <= t <= this.shift+this.a)
                return evalpoly(t-this.shift, this.n)
            elseif (this.shift-this.a <= t < this.shift)
                if this.sym == 'e'
                    return conj(evalpoly(this.shift-t, this.n))
                elseif this.sym == 'o'
                    return -1*conj(evalpoly(this.shift-t, this.n))
                end
            else
                return 0
            end
        end
    end
    if typeof(this) == collection
        result = 0
        for el in this.elements
            result+=evaluate(el, t)
        end
        return result
    end
    return 0
end
# add function evalDer to get derivative N of collection at a point
#  from left or from right with an (optional?) argument
#  if not specified, choose a direction, print error if left/right do not match


function evalSpectrum(this::Union{collection,element}, w)
    if typeof(this) == element
        if w == 0
            result_left = big(0)
            coeff = big(this.a)
            for i in 0:this.N
                result_left += coeff*this.u[i+1]
                coeff*=this.a//(i+2)
            end
            for i in (this.N+1):(2*this.N+1)
                result_left += coeff*this.h[i-this.N]
                coeff*=this.a//(i+2)
            end
            result_right = big(0)
            coeff = big(this.a)
            for i in 0:this.N
                result_right += coeff*this.v[i+1]
                coeff*=-1*this.a//(i+2)
            end
            for i in (this.N+1):(2*this.N+1)
                result_right += coeff*this.k[i-this.N]
                coeff*=-1*this.a//(i+2)
            end
            
            if result_left != result_right
                print("error: result_left = ", result_left, "\n       result_right = ", result_right, "\n")
            end
            if this.sym == 'e'
                result_left += conj(result_left)
            elseif this.sym == 'o'
                result_left -= conj(result_left)
            end
            return result_left
        else
            # effects of symmetry are already computed in this.spec
            result = 0
            for term in this.spec.terms
                result += (term.coeff+term.eiwCoeff*exp(big(-1.0im)*term.a*w)+term.cosCoeff*cos(big(this.a)*w)+term.sinCoeff*sin(big(this.a)*w))/big(w)^term.wDeg
            end
            result *= exp(big(-1.0im)*this.shift*w)
            return result
        end
    end
    if typeof(this) == collection
        if isempty(this.spec.terms)
            constructSpectrum(this)
        end
        if w == 0
            # the spectrum as w->0 can only be computed on elements
            result = 0
            for el in this.elements
                result += evalSpectrum(el, w)
            end
            return result
        else
            result = 0
            for term in this.spec.terms
                result += (term.coeff+term.eiwCoeff*exp(big(-1.0im)*term.a*w)+term.cosCoeff*cos(big(term.a)*w)+term.sinCoeff*sin(big(term.a)*w))*exp(big(-1.0im)*term.shift*w)/big(w)^term.wDeg
            end
            return result
        end
    end
end

using Plots

# pass by reference recursive function to find lowest and highest values of t
# [lowt] and [hight] are how the variables are passed in
function getBounds(this::Union{collection,element}, lowt, hight)
    if typeof(this) ==  collection
        for el in this.elements
            getBounds(el, lowt, hight)  # variables are already references
        end
    elseif typeof(this) == element
        lowthis = (this.sym == 'n' ? this.shift : this.shift - this.a)
        if lowthis < lowt[1]
            lowt[1] = lowthis
        end
        highthis = this.shift + this.a
        if highthis > hight[1]
            hight[1] = highthis
        end
    end
end
            
function plotElement(this::Union{collection,element})
    # hardcoded points per unit interval until later
    ppui = 50
    lowt = [0]
    hight = [0]
    getBounds(this, lowt, hight)
    t = lowt[1]:(1.0/ppui):hight[1]
    fcn = zeros(Float64, length(t))
    for i in 1:length(t)
        fcn[i] = evaluate(this, t[i])
    end
    plot(t, fcn)
end

function plotSpectrum(this::Union{collection,element})
    # hardcoded num_decades and num_points to plot
    num_decades = 2
    num_points = 400  
    lowf = (0.5*10^(-1*num_decades/2))
    highf = (0.5*10^(num_decades/2))
    inc = (highf - lowf) / (num_points - 1)
    f = lowf:inc:highf
    omega = 2*pi*f
    fcn = zeros(Float64, length(omega))
    highest = 0
    for i in 1:length(omega)
        fcn[i] = abs(evalSpectrum(this, omega[i]))
        if fcn[i] < highest
            if fcn[i] <= (10.0^(-6)*highest)
                fcn[i] = 10.0^(-6)*highest
            end
        else
            highest = fcn[i]
        end
    end
    plot(f, fcn, xaxis = :log, yaxis = :log)
end

# Next steps:
#  Add scalar multiply function on elements/collections
#  Create solver for the critical problems
#    --learn if/how to join vectors and matrices 
#  Apply to truncated sinc 
#    --needs method to numerically integrate truncated sinc for lim w->0 F(w)
#    --needs method to calculate Nth derivative of sinc at truncation points
#  Make better plotting 
#
#  math ideas:
#  multiplication (has matrix with binomial coefficients)
#    --needs partitioning scheme
#    --will need function to join terms with same N, shift, a
#    --may also need embedding function listed below
#  convolution (like mulitplication but adds integration)
#  inner product (multiplication and integration)
#  orthogonalization
#  embedding, function that raises N to higher number
