export getBounds, plotElement, plotSpectrum

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
