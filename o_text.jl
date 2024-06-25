export complexToString, strSpectrum, vectorToString, strElement

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
