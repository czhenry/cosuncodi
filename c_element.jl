# the spectrum struct stores information about the Fourier Transform
# of an element or a collection of elements

export spectralTerm, nullTerm, spectrum, nullSpectrum, element, collection
export nullCollection, sumSpectrums, reduceSpectrum, constructSpectrum

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

# the element struct stores information about a polynomial of degree at most 2N+1
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
    N=N%4
    if N<0
        N=N+4
    end
    N<2 ? (N==0 ? 1 : im) : (N==2 ? -1 : -im)
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
