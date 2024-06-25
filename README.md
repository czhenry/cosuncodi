# cosuncodi
CoSuNCoDi v0.0.1 

Compactly Supported N-times Continuously Differentiable Piecewise Polynomial Basis using 

the Julia programming language


usage:

<code> julia -i cosuncodi.jl </code>

or in scripts

<code> include("cosuncodi.jl") </code>


This library includes:

1.  data structures:

element

collection

spectrum


2.  element constructors:

u

v

w


3.  basis collection constructors

uBasis

vBasis


and 4. functions to plot collections/elements and spectrums and to print strings in LaTeX

plotSpectrum

plotElement

strSpectrum

StrElement




using which one can represent piecewise polynomials by specifying their value and the value

of their first N derivatives at each of the left and right end points of an interval.


The element and collection data structures each have a spectrum that represents the Fourier

Transform of it.



