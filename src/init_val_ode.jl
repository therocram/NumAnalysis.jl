##########################################################
# Initial Value ODE
# File of functions used for approximating the solutions to
# well-posed initial value ODE problems
#
# Author: Masen Pitts
# Last Modified: 2/16/2022 (MM/DD/YYYY)
# Version: 1.0
#########################################################

# Format the output table initial value problem approximation with the initial condition 
# y(t0) = w0
# using the whitespace tokens "token1" and "token2"
#
# The true solution "y" is optionally passed to the function to compute and print
# absolute error
function initprint(t0, w0, token1, token2, y=nothing)
    if y!== nothing
        println(string("i", token2, "t_i", token1, "w_i", token1, "|y(t_i) - w_i|"))
    else
        println(string("i", token2, "t_i", token1, "w_i"))
    end

    printstep(0, t0, w0, token2, y)
    
    # Print inital conditions
end

# Output the "i"th step of an initial value problem approximation with time "t"
# and solution value "w", using the whitespace tokens "token1" and "token2"
#
# The true solution "y" is optionally passed to the function to compute and print
# absolute error
function printstep(i, t, w, token2, y=nothing)
    if y !== nothing
        println(string(i, token2, @sprintf("%.6f", t), token2, @sprintf("%.6f", w), 
                        token2, @sprintf("%.6f", abs(y(t) - w))))
    else
        println(string(i, token2, @sprintf("%.6f", t), token2, @sprintf("%.6f", w)))
    end
end

# Initializes an IVP for various approximation methods given the interval 
#   a <= t <= b
# the initial condition y(a) = y0, and the total number of steps N
#
# The true solution "y" is optionally given if the error of the approximation
# is to be printed
function ivpsetup(a, b, y0, N, y=nothing, printres=false)
    # Determine spacing of mesh points
    h = (b-a)/N

    t = zeros(N+1)
    w = zeros(N+1)

    # Set up initial values
    t[1] = a
    w[1] = y0
    
    # Special whitespace strings used for table formatting
    token1 = "        "
    token2 = "   "

    # Format output table and print inital conditions
    if printres
        initprint(t[1], w[1], token1, token2, y)
    end

    return h, t, w, token2
end

# Add this to documentation soon
# See also [`rungekuttaO4`](@ref), [`euler`](@ref), [`midpoint`](@ref)

"""
    ivpsolve(f, a, b, y0, N; solver=rungekuttaO4, y=nothing, printres=false)

Approximate the solution ``y(t)`` to the initial-value problem
```math
y' = f(t,y), \\hspace{0.2cm} a \\leq t \\leq b, \\hspace{0.2cm} y(a) = y_0
```
at `N` + 1 equally spaced values of ``t`` on the interval [`a`,`b`] using the
one step solver method given by `solver`.

The function returns a tuple of the form `(t, w)`, where `t` and `w` are Vectors containing all of the 
independent and dependent variable values in the solution set, respectively. By default, 
the Runge-Kutta method of order 4 (RKO4) is used:

```jldoctest
julia> f(t,y) = -2t*y
f (generic function with 1 method)

julia> ivpsolve(f, 0, 1, 1, 5)
([0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0], [1.0, 0.9607893333333333, 0.8521429680674133, 0.6976755803411453, 0.527297771054652, 0.3679036697909508])
```

The user can optionally have the method print the results of each step of the solver by 
setting `printres = true`. Note that in the example below the actual function output is suppressed
by the `;`:

```jldoctest
julia> ivpsolve(f, 0, 1, 1, 5, printres=true);
i   t_i        w_i
0   0.000000   1.000000
1   0.200000   0.960789
2   0.400000   0.852143
3   0.600000   0.697676
4   0.800000   0.527298
5   1.000000   0.367904
```

Providing the true solution `y` will automatically print the results of each step along
with the true error of the approximation at each step:

```jldoctest
julia> y(t) = exp(-t^2)
y (generic function with 1 method)

julia> ivpsolve(f, 0, 1, 1, 5, y=y);
i   t_i        w_i        |y(t_i) - w_i|
0   0.000000   1.000000   0.000000
1   0.200000   0.960789   0.000000
2   0.400000   0.852143   0.000001
3   0.600000   0.697676   0.000001
4   0.800000   0.527298   0.000005
5   1.000000   0.367904   0.000024
```

In addition to RKO4 there are two other single-step solver methods (Euler's method and 
the Midpoint/Runge-Kutta Order 2 method). Out of the three, ROK4 will generally be the most
accurate. Using Euler's Method, for example, will almost always give far less accurate results:

```jldoctest
julia> ivpsolve(f, 0, 1, 1, 5, solver=euler, y=y);
i   t_i        w_i        |y(t_i) - w_i|
0   0.000000   1.000000   0.000000
1   0.200000   1.000000   0.039211
2   0.400000   0.920000   0.067856
3   0.600000   0.772800   0.075124
4   0.800000   0.587328   0.060036
5   1.000000   0.399383   0.031504
```
"""
function ivpsolve(f, a, b, y0, N; solver=rungekuttaO4, y=nothing, printres=false)
    # If y(t), the true solution, is provided then the results will automatically
    # be printed.
    if printres == false && y !== nothing
        printres = true
    end
    
    # Setup mesh spacing, solution arrays, and custom whitespace token
    h, t, w, token2 = ivpsetup(a, b, y0, N, y, printres)

    # Approximate y(t) at each mesh point
    for i in 2:N+1
        t[i], w[i] = solver(t[i-1], w[i-1], f, h)

        if printres
            printstep(i-1, t[i], w[i], token2, y)
        end
    end

    return t, w
end

# Euler's Method
function euler(t, y, f, h)
    tnext = t + h 
    ynext = y + h*f(t,y)

    return tnext, ynext  
end

# Runge-Kutta (Order Four)
function rungekuttaO4(t, y, f, h)
    k1 = h*f(t,y)
    k2 = h*f(t+0.5*h, y+0.5*k1)
    k3 = h*f(t+0.5*h, y+0.5*k2)
    k4 = h*f(t+h, y+k3)

    tnext = t + h
    ynext = y + (k1 + 2*k2 + 2*k3 + k4)/6

    return tnext, ynext
end

# Midpoint method (Runge-Kutta Order 2)
function midpoint(t, y, f, h)
    k = y + 0.5*h*f(t,y)

    tnext = t + h
    ynext = y + h*f(t+0.5*h,k)

    return tnext, ynext
end

# # Euler's Method
# # Approximate the solution y(t) to the initial-value problem
# #   y' = f(t,y), a <= t <= b, y(a) = y0
# # at N + 1 equally spaced numbers on the interval [a,b]
# #
# # The true solution "y" is optionally given if the error of the approximation
# # is to be printed
# function euler(f, a, b, y0, N, y=nothing)
#     h, t, w, token2 = ivpsetup(a, b, y0, N, y)

#     # Approximate y(t) at each mesh point
#     for i in 2:N+1
#         w[i] = w[i-1] + h*f(t[i-1],w[i-1])
#         t[i] = t[i-1] + h

#         printstep(i-1, t[i], w[i], token2, y)
#     end

#     return t, w
# end

# # Runge-Kutta (Order Four)
# # Approximate the solution y(t) to the initial-value problem
# #   y' = f(t,y), a <= t <= b, y(a) = y0
# # at N + 1 equally spaced numbers on the interval [a,b]
# #
# # The true solution "y" is optionally given if the error of the approximation
# # is to be printed
# function rungekuttaO4(f, a, b, y0, N, y=nothing)
#     h, t, w, token2 = ivpsetup(a, b, y0, N, y)

#     # Approximate y(t) at each mesh point
#     for i in 2:N+1
#         t1 = t[i-1]
#         w1 = w[i-1]

#         k1 = f(t1,w1)
#         k2 = f(t1+0.5*h, w1+0.5*k1)
#         k3 = f(t1+0.5*h, w1+0.5*k2)
#         k4 = f(t1+h, w1+k3)

#         w[i] = w1 + h*(k1+k2+k3+k4)/6
#         t[i] = t1 + h

#         printstep(i-1, t[i], w[i], token2, y)
#     end

#     return t, w
# end
