#########################################################
# Calculus
# File of functions numerically approximating derivatives
# and integrals

# Author: Masen Pitts
# Last Modified: 4/2/2022 (MM/DD/YYYY)
#########################################################

#using FastGaussQuadrature

# Composite Simpson's Rule.
# Approximates the integral I = Int(f(x), a < x < b) using
# Simpson's Rule on "n" subintervals.
#
# NOTE: Using this function with n=2 is equivalent to performing
# a single Simpson's Rule approximation over [a,b]
function simpson(f, a, b, n)
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


"""
    adaptivequad(f, a, b, tol, N, base=true)

Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

within a given tolerance `tol` using a recursive adaptive quadrature 
method based on Simpson's Rule for a maximum of `N` subintervals.

The parameter `base` is a boolean value that controls when the 
program initializes and when it fails. It should always be set to `true`
when being called by the user.

# Examples


See also [`simpson`](@ref), [`gaussianquad`](@ref), [`trapezoid`](@ref)
"""
function adaptivequad(f, a, b, tol, N, base=true)
    # Initializes subinterval tracking variable "n" at base iteration
    if base 
        global n = 1 
    end

    # Fail condition
    if n > N
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
                println("\nI = ", left + right, "\nSubintervals required: ", n)
                return left + right
            else
                return left + right
            end
        end
    end
end


"""
    gaussianquad(f::Function, a, b, n; sub = 1)

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
is 0.804776489.
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

See also [`adaptivequad`](@ref), [`simpson`](@ref), [`trapezoid`](@ref)
"""
function gaussianquad(f::Function, a, b, n; sub=1)
    # Use FastGaussQuadrature package to very quickly obtain nodes and weights
    # of Legendre polynomial expansion
    nodes, weights = FastGaussQuadrature.gausslegendre(n)

    # Divide region of integration into "sub" number of subintervals
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
            subsum += weights[i]*f( ( (b-a)*nodes[i] + (b+a) )/2 )
        end

        sum += subsum*(b-a)/2
    end

    return sum
end

# This code is obsolete compared to the FastGaussQuadrature package's
# capabilites. Still a cool learning exercise. 
#
# using Polynomials, SpecialPolynomials

# """
#     getlegendre(n)

# Return the nth order Legendre polynomial P_n(x) as implemented
# by the `SpecialPolynomials` packages.
# """
# function getlegendre(n)
#     last = zeros(n+1); last[end] = 1

#     return Legendre(last)
# end

# """
#     getweight(nodes, i)

# Return the Gaussian quadrature weight corresponding to the ith node x_i, 
# where the P(x_i) = 0 and P(x) is the appropriate Legendre polynomial.
# """
# function getweight(nodes, i)
#     xi = nodes[i]

#     p = Polynomial(1)

#     for j in eachindex(nodes)
#         if j != i
#             nodediff = nodes[i] - nodes[j]
#             p *= Polynomial([-nodes[j]/(nodediff), 1/(nodediff)])
#         end
#     end
    
#     return integrate(p, -1, 1)
# end

# """
#     gausslegendre(n)

# Return the nodes and weights corresponding to a
# Gaussian Quadrature approximation of order n in the Legendre polynomials.

# Values of 1 ≤ n ≤ 20 should generally be used. Larger values of n
# can be used, but round off error begins to become non-negligible.

# # Examples
# ```jldoctest
# julia> nodes, weights = gausslegendre(4);

# julia> nodes
# 4-element Vector{Float64}:
#  -0.8611363115940526
#  -0.33998104358485626
#   0.33998104358485626
#   0.8611363115940526

# julia> weights
# 4-element Vector{Float64}:
#  0.6521451548625463
#  0.3478548451374537
#  0.6521451548625463
#  0.3478548451374537
# ```
# """
# function gausslegendre(n)
#     P = getlegendre(n)

#     nodes = real(roots(P)) # Most computationally demanding part
#     weights = zeros(n)

#     if n % 2 == 0
#         half = div(n,2) + 1:n
#         halfweights = [getweight(nodes, i) for i in half]
#         weights[half] = halfweights
#         weights[1:div(n,2)] = reverse(halfweights)
#     else
#         midhalf = div(n,2) + 1:n
#         halfweights = [getweight(nodes, i) for i in midhalf]
#         weights[midhalf] = halfweights
#         weights[1:div(n,2)] = reverse(halfweights)[1:end-1]
#     end

#     return nodes, weights
# end
