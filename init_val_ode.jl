##########################################################
# Initial Value ODE
# File of functions used for approximating the solutions to
# well-posed initial value ODE problems
#
# Author: Masen Pitts
# Last Modified: 2/16/2022 (MM/DD/YYYY)
# Version: 1.0
#########################################################

using Printf

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

# Euler's Method
# Approximate the solution y(t) to the initial-value problem
#   y' = f(t,y), a <= t <= b, y(a) = y0
# at N + 1 equally spaced numbers on the interval [a,b]
#
# The true solution "y" is optionally given if the error of the approximation
# is to be printed
function euler(f, a, b, y0, N, y=nothing)
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
    initprint(t[1], w[1], token1, token2, y)

    # Approximate y(t) at each mesh point
    for i in 2:N+1
        w[i] = w[i-1] + h*f(t[i-1],w[i-1])
        t[i] = t[i-1] + h

        printstep(i-1, t[i], w[i], token2, y)
    end

    return t, w
end