##########################################################
# Calculus
# File of helper functions numerically approximating derivatives
# and integrals

# Author: Masen Pitts
# Last Modified: 11/05/2022 (MM/DD/YYYY)
# Version: 1.1
#########################################################


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

# Adaptive Quadrature Method
# Approximates the integral I = Int(f(x), a < x < b) within
# a given tolerance "tol" using a recursive adaptive quadrature method
# based on Simpson's Rule with n=4 for a maximum of N subintervals
#
# The parameter "base" is a boolean value that controls when the 
# program initializes and when it fails. It should always be set to "true"
# when being called by the user.
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
