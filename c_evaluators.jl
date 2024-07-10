export evaluate, evalSpectrum
export evalDer

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

function doNthDer(this::element, t, N)
    if N < 0
        print("error: evalDer does not do derivatives with N<0\n")
        return 0
    elseif N == 0
        return evalpoly(t-this.shift, this.n)
    else
        derPoly = deepcopy(this.n)
        # note to avoid read after write errors order from left to right
        for i in 1:N
            for j in 2:(2*this.N+2)
                derPoly[j-1]=(j-1)*derPoly[j]
            end
            derPoly[2*this.N+2] = 0
        end
        return evalpoly(t-this.shift, derPoly)
    end
end

function evalDer(this::Union{collection,element}, t, N=1, dir='r')
    if typeof(this) == element
        if this.sym == 'n'
            if (this.shift < t < this.shift+this.a)
                # we're inside an open set, direction doesn't matter
                return doNthDer(this, t, N)
            elseif t == this.shift
                return (dir == 'r' ) ? doNthDer(this, t, N) : 0
            elseif t == this.shift+this.a
                return (dir == 'l' ) ? doNthDer(this, t, N) : 0
            else
                return 0
            end
        else
            # We must be symmetric
            if (this.shift < t < this.shift+this.a)
                return doNthDer(this, t, N)
            elseif (this.shift-this.a < t < this.shift)
                if this.sym == 'e'
                    return (iseven(N) ? 1 : -1)*conj(doNthDer(this, t, N))
                elseif this.sym == 'o'
                    return (isodd(N) ? 1 : -1)*conj(doNthDer(this, t, N))
                end
            elseif t == this.shift+this.a
                return (dir == 'l' ) ? doNthDer(this, t, N) : 0
            elseif t == this.shift-this.a
                if dir == 'r'
                    if this.sym == 'e'
                        return (iseven(N) ? 1 : -1)*conj(doNthDer(this, t, N))
                    elseif this.sym == 'o'
                        return (isodd(N) ? 1 : -1)*conj(doNthDer(this, t, N))
                    end
                else
                    return 0
                end
            elseif t == this.shift
                if dir == 'r'
                    return doNthDer(this, t, N)
                else
                    if this.sym == 'e'
                        return (iseven(N) ? 1 : -1)*conj(doNthDer(this, t, N))
                    elseif this.sym == 'o'
                        return (isodd(N) ? 1 : -1)*conj(doNthDer(this, t, N))
                    end
                end
            else
                return 0
            end
        end
    end
    if typeof(this) == collection
        result = 0
        for el in this.elements
            result+=evalDer(el, t, N, dir)
        end
        return result
    end
    return 0
end



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
