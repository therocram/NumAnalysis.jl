#########################################################
# Calculus
# File of functions numerically approximating derivatives
# and integrals

# Author: Masen Pitts
# Last Modified: 4/2/2022 (MM/DD/YYYY)
#########################################################

"""
    simpson(f::Function, a, b, n=2)

Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

using a composite Simpson's Rule over `n` subintervals. Using this
method with `n = 2` is equivalent to performing a single Simpson's
Rule approximation over [`a`,`b`].

# Examples
The true value of the integral in the example below accurate to 16
decimal places is 0.3386826421694169:
```jldoctest simpson
julia> f(x) = x^3 * cosh(x);

julia> simpson(f, 0, 1, 2)
0.35114893623640564
```

In general, increasing `n` increases the accuracy of the approximation
as long as `n` is even:
```jldoctest simpson
julia> simpson(f, 0, 1, 10)
0.33870448867643016
```

The method will not give accurate results for odd `n` and will warn the
user automatically if the inputted `n` is odd:
```jldoctest simpson
julia> simpson(f, 0, 1, 7)
┌ Warning: Method will not give accurate results for odd n
└ @ NumAnalysis ~/NumAnalysis.jl/src/calculus.jl:54
0.2847753864365918
```

See also [`adaptivequad`](@ref), [`trapezoid`](@ref), [`gaussianquad`](@ref)
"""
function simpson(f::Function, a, b, n=2)
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

"""
    trapezoid(f::Function, a, b, n=1)

Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

using a composite Trapezoidal Rule over `n` subintervals. Using this
method with `n = 1` is equivalent to performing a single Trapezoidal
Rule approximation over [`a`,`b`].

# Examples
The true value of the integral in the example below accurate to 16
decimal places is 0.3386826421694169:
```jldoctest trapezoid
julia> f(x) = x^3 * cosh(x);

julia> trapezoid(f, 0, 1, 1)
0.7715403174076219
```

In general, increasing `n` increases the accuracy of the approximation:
```jldoctest trapezoid
julia> trapezoid(f, 0, 1, 10)
0.34351419965225494
```

See also [`simpson`](@ref), [`gaussianquad`](@ref)
"""
function trapezoid(f::Function, a, b, n=1)
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
method applying Simpson's Rule (see [`simpson`](@ref)) over a maximum of `N` subintervals.

The parameter `base` is a boolean value that tracks the initial call 
of the method and controls when it fails. It should always be set to `true`
when being called by the user.

# Examples
The method can obtain reasonably accurate results with very few subintervals.
The true value of the integral in the example below accurate to 20 decimal
places is 0.16060279414278839202:
```jldoctest adaptivequad
julia> f(x) = x^2 * exp(-x);

julia> adaptivequad(f, 0, 1, 1e-5, 10)

Subintervals required: 3
0.16060529683648378
```
However, the number of required calculations grows quickly when very high
precision is required:
```jldoctest adaptivequad
julia> adaptivequad(f, 0, 1, 1e-17, 5000)

Subintervals required: 3186
0.16060279414278839
```

See also [`simpson`](@ref), [`adaptivegaussquad`](@ref)
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
        if base println("\nSubintervals required: ", n) end
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
    gaussianquad(f, a, b, n; sub=1, nwlist=nothing, checklist=true)
    
Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

using a `n` order Gaussian quadrature approximation on `sub`
subintervals. The quadrature nodes and weights are optionally
provided as a size 2 tuple of Vectors or Vector of Vectors with
`nwlist`. 

By default, the method will check `nwlist` to make sure
it has appropriate type and size. These safety checks can be skipped
by setting `checklist = false`.

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
```jldoctest gaussianquad
julia> f(x) = sin(x^2);

julia> gaussianquad(f, 0, 2, 5)
0.8048689412592385
```
Accuracy can be further improved by making use of subintervals:
```jldoctest gaussianquad
julia> gaussianquad(f, 0, 2, 4, sub=4)
0.8047764537878016
```

See also [`adaptivequad`](@ref), [`adaptivegaussquad`](@ref)
"""
function gaussianquad(f, a, b, n; sub=1, nwlist=nothing, checklist=true)
    # Use FastGaussQuadrature package to very quickly obtain nodes and weights
    # of Legendre polynomial expansion
    if nwlist === nothing
        nodes, weights = FastGaussQuadrature.gausslegendre(n)
    # Allows optional passing of nodes and weights for more efficient repeated calls
    elseif !checklist
        nodes, weights = nwlist
    # Ensures that the passed list of quadrature nodes and weights is either a Tuple of 2 Vectors
    # or a Vector of 2 Vectors. This check is skipped if checklist=false
    elseif isa(nwlist, Tuple{Vector,Vector}) || 
                ( isa(nwlist, Vector{Vector{Float64}}) && size(nwlist)[1] == 2 )
        
        nodes, weights = nwlist

        # There must be an equal number of quadrature nodes and weights
        if size(nodes) != size(weights)
            return @error("The passed lists of quadrature nodes and weights are of different sizes", nwlist)
        # The number of nodes and weights must match the order of the approximation
        elseif size(nodes)[1] != n
            return @error("The passed lists of quadrature nodes and weights must be of size n", n, nwlist)
        end
    else
        return @error("The passed list of quadrature nodes and weights must be of type Tuple{Vector,Vector} or 
        a Vector{Vector{Float64}} with a size of 2", nwlist)
    end

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

"""
    adaptivegaussquad(f, a, b, n, tol, submax, base=true)

Approximates the integral 

```math
\\int_{a}^{b} f(x) dx 
```

within a given tolerance `tol` using a recursive adaptive quadrature 
method applying order `n` Gaussian Quadrature approximations (see [`gaussianquad`](@ref)) 
over a maximum of `submax` subintervals.

The parameter `base` is a boolean value that tracks the initial call 
of the method and controls when it fails. It should always be set to `true`
when being called by the user.

# Examples
Even for small values of `n` this method can obtain highly accurate results with
fewer subintervals than the Simpson's Rule analogue (see [`adaptivequad`](@ref)).
The true value of the integral in the example below accurate to 20 decimal
places is 0.16060279414278839202:
```jldoctest adapgauss
julia> f(x) = x^2 * exp(-x);

julia> adaptivegaussquad(f, 0, 1, 2, 1e-17, 3000)

Subintervals required: 2574
0.16060279414278839
```
The number of required subintervals drastically decreases for even slightly larger
values of `n`:
```jldoctest adapgauss
julia> adaptivegaussquad(f, 0, 1, 4, 1e-17, 10)

Subintervals required: 9
0.16060279414278839
```

See also [`gaussianquad`](@ref), [`adaptivequad`](@ref)
"""
function adaptivegaussquad(f, a, b, n, tol, submax, base=true)
    # Initializes crucial tracking variables and constants at base iteration
    if base 
        # Subinterval tracking variable
        global subnum = 1
        # Constant used to approximate the error of each quadrature approximation
        global errormod = floor( (2/3)*(2^(2*n) - 1) )
        # List of untransformed gaussian quadrature nodes and weights
        global NWLIST = 
            FastGaussQuadrature.gausslegendre(n)::Tuple{Vector{Float64},Vector{Float64}}
    end

    # Fail condition
    if subnum::Int64 > submax
        return nothing
    end

    # Approximates I with sub=1 and sub=2 nth order Gaussian Quadrature calculations
    G = gaussianquad(f, a, b, n, nwlist=NWLIST, checklist=false)
    G2 = gaussianquad(f, a, b, n, sub=2, nwlist=NWLIST, checklist=false)

    # Recursively apply method to left and right subintervals when 
    # approximations do not meet a defined tolerance
    if abs(G2 - G) < errormod*tol
        if base println("\nSubintervals required: ", subnum) end
        return G2
    else
        global subnum += 1
        mid = (a + b)/2

        left = adaptivegaussquad(f, a, mid, n, tol/2, submax, false)
        right = adaptivegaussquad(f, mid, b, n, tol/2, submax, false)

        # Neatly fail the method if the number of subintervals exceeds
        # the allowed maximum and return a proper value otherwise. 
        if left === nothing || right === nothing
            if base
                return @error("Method failed: exceeded max subintervals", submax)
            else
                return nothing
            end
        else
            if base
                println("\nSubintervals required: ", subnum)
                return left + right
            else
                return left + right
            end
        end
    end
end
