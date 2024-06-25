export evaluate, evalSpectrum

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
