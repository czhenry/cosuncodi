# cosuncodi
CoSuNCoDi v0.0.1 

Compactly Supported N-times Continuously Differentiable Piecewise Polynomial Basis using 

the Julia programming language


usage:

<code> julia -i element.jl </code>

or in scripts

<code> include("element.jl") </code>


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


This first version of the README.md file and code are expected to change.  This is just to 

note the first *working* version of the code, before re-organizing it.

