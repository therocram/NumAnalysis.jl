#########################################################
# Calculus
# File of functions numerically approximating derivatives
# and integrals

# Author: Masen Pitts
# Last Modified: 4/2/2022 (MM/DD/YYYY)
#########################################################

# Composite Simpson's Rule.
# Approximates the integral I = Int(f(x), a < x < b) using
# Simpson's Rule on "n" subintervals.
#
# NOTE: Using this function with n=2 is equivalent to performing
# a single Simpson's Rule approximation over [a,b]
function simpson(f::Function, a, b, n)
    # Give warning if n is not even
    if n % 2 != 0
        @warn("Method will not give accurate results for odd n")
    end

    h = (b - a)/n

    y0 = f(a) + f(b)
    y1 = 0; y2 = 0

    for i in 1:n-1
        y = a + i*h

        if i % 2 == 0
            y2 += f(y)
        else
            y1 += f(y)
        end
    end

    return (h/3)*(y0 + 2*y2 + 4*y1)
end

# Composite Trapezoidal Rule.
# Approximates the integral I = Int(f(x), a < x < b) using
# the Trapezoidal Rule on "n" subintervals.
#
# NOTE: Using this function with n=1 is equivalent to performing
# a single Trapezoidal Rule approximation over [a,b]
function trapezoid(f, a, b, n)
    h = (b - a)/n

    y0 = f(a) + f(b)
    y1 = 0

    for i in 1:n-1
        y = a + i*h
        y1 += f(y)
    end

    return h*(y0 + 2*y1)/2
end

# Add these to both docs later
#, [`simpson`](@ref), [`trapezoid`](@ref)


"""
    adaptivequad(f, a, b, tol, N, base=true)

Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

within a given tolerance `tol` using a recursive adaptive quadrature 
method applying Simpson's Rule over a maximum of `N` subintervals.

The parameter `base` is a boolean value that tracks the initial call 
of the method and controls when it fails. It should always be set to `true`
when being called by the user.

# Examples
The method can obtain reasonably accurate results with very few subintervals.
The true value of the integral in the example below accurate to 20 decimal
places is 0.16060279414278839202:
```jldoctest
julia> f(x) = x^2 * exp(-x)
f (generic function with 1 method)

julia> adaptivequad(f, 0, 1, 1e-5, 10)

Subintervals required: 3
0.16060529683648378
```
However, the number of required calculations grows quickly when very high
precision is required:
```jldoctest
julia> adaptivequad(f, 0, 1, 1e-17, 5000)

Subintervals required: 3186
0.16060279414278839
```

See also [`gaussianquad`](@ref)
"""
function adaptivequad(f, a, b, tol, N, base=true)
    # Initializes subinterval tracking variable "n" at base iteration
    if base 
        global n = 1 
    end

    # Fail condition
    if n::Int64 > N
        return nothing
    end

    # Approximates I with n=2 and n=4 Simpson's Rule calculations
    S = simpson(f, a, b, 2)
    S4 = simpson(f, a, b, 4)

    # Recursively apply method to left and right subintervals when 
    # approximations do not meet a defined tolerance
    if abs(S - S4) < 10*tol
        if base println("\nI = ", S4, "\nSubintervals required: ", n) end
        return S4
    else
        global n += 1
        mid = (a + b)/2

        left = adaptivequad(f, a, mid, tol/2, N, false)
        right = adaptivequad(f, mid, b, tol/2, N, false)

        # Neatly fail the method if the number of subintervals exceeds
        # the allowed maximum and return a proper value otherwise. 
        # If you want the method to fail AS SOON
        # as any of the recursive iterations fail, replace the code with
        # the logic below (this will give a much messier error message):
        # if left === nothing || right === nothing
        #     error("\nMethod failed: exceeded N subintervals. N = ", N)
        # else
        #     return left + right
        if left === nothing || right === nothing
            if base
                return @error("Method failed: exceeded N subintervals", N)
            else
                return nothing
            end
        else
            if base
                println("\nSubintervals required: ", n)
                return left + right
            else
                return left + right
            end
        end
    end
end


"""
    gaussianquad(f, a, b, n; sub = 1)

Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

using a `n` order Gaussian quadrature approximation on `sub`
subintervals.

By default, the Gaussian quadrature nodes are spread over the
given interval (`a`,`b`), meaning that the error of the approximation
is reduced by using a higher `n`.

However, error can also be significantly reduced by applying the same
order of Gaussian quadrature to a specified number of subsets of the
original interval (`a`,`b`). In the case where `sub` ``\\neq`` 1, 
Gaussian quadrature is seperately applied to `sub` subintervals, each
of width `h = (b - a)/sub`.

# Examples
Even with a very low order `n` we can still produce accurate results.
The true value of the integral estimated below to 9 decimal places 
is 0.804776489:
```jldoctest
julia> f(x) = sin(x^2)
f (generic function with 1 method)

julia> gaussianquad(f, 0, 2, 5)
0.8048689412592385
```
Accuracy can be further improved by making use of subintervals:
```jldoctest
julia> gaussianquad(f, 0, 2, 4, sub=4)
0.8047764537878016
```

See also [`adaptivequad`](@ref)
"""
function gaussianquad(f, a, b, n; sub=1)
    # Use FastGaussQuadrature package to very quickly obtain nodes and weights
    # of Legendre polynomial expansion
    nodes, weights = FastGaussQuadrature.gausslegendre(n)

    # Divide region of integration into "sub" number of equally spaced
    # subintervals
    h = (b-a)/sub
    leftends = [a + j*h for j in 0:sub-1] 
    rightends = [a + j*h for j in 1:sub]

    sum = 0 # Accumulates integral area

    # Estimate integral over each subinterval, then return the total
    # sum
    for j in 1:sub
        a = leftends[j]
        b = rightends[j]

        # Perfrom Gaussian Quadrature integration by transforming
        # integration variable to a region on the interval [-1,1]
        subsum = 0
        for i in eachindex(nodes)
            subsum += weights[i]*f( ( h*nodes[i] + (b+a) )/2 )
        end

        sum += subsum*h/2
    end

    return sum
end
