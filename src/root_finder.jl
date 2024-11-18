#########################################################
# Root Finder
# File of functions used for finding the roots and
# fixed points of real functions of one variable

# Author: Masen Pitts
# Last Modified: 11/05/2022 (MM/DD/YYYY)
# Version: 1.1
#########################################################

# Bisection Method
# Implements the bisection method on the interval [a,b]
# to approximate the root of a function "f" within some tolerance
# "tol" for a maximum of "N" iterations.
"""
    bisecmethod(f::Function, a, b, tol, N)

Uses the Bisection Method to approximate the root of a single-variable 
function `f` on the interval [`a`,`b`] within the set tolerance `tol`
for a maximum of `N` iterations.
"""
function bisecmethod(f::Function, a, b, tol, N)
    # Checks for roots at the endpoints of [a,b]
    if f(a) == 0
        return a
    elseif f(b) == 0
        return b
    end

    # Prints header of output table
    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
    "Error")
    println("--------------------------------------------\
     ---------------")

    # Throws error if f has no root on [a,b]
    if f(a)f(b) > 0
        return @error("function values at endpoints must \
                have different signs", f(a), f(b))
    end

    # Keeps track of previous root approximation value
    lastp = a

    for i = 1:N
        # New guess for root value is the midpoint of the current interval
        p = a + (b-a)/2
        val = f(p)

        # Checks if new root guess is actually a root
        if val == 0
            println("\nSolution at ", p, ". f(", p, ") = 0")
            return p
        end

        # Estimates error using difference between successive root guesses
        error = p == 0 ? abs(p-lastp) : abs(p-lastp)/abs(p) 

        # Print results for current iteration
        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, val, error)

        # Return current guess if estimated error is within tol
        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end

        # Selects the half of the previous interval that still has a root of f
        # for the next iteration
        if sign(f(a))*sign(val) > 0
            a = p
        else
            b = p
        end

        # Update value of previous root guess
        lastp = p
    end

    return @error("method failed after N interations. N = ", N)
end

# Fixed Point
# Approximates the fixed point of some function "g" given
# an inital guess "p0" within a tolerance of "tol" for a 
# maximum of "N" iterations.
function fixedpoint(g, p0, tol, N)
    for i in 1:N
        # New guess for root is g evaluated at the previous guess
        p = g(p0)
        # Print result once estimated error falls within tol
        if p == 0 ? abs(p-p0) < tol : abs(p-p0)/abs(p) < tol
            println("Required Iterations: ", i)
            return p
        end

        # Update value of previous root guess
        p0 = p
    end

    return @error("method failed after N interations", N)
end

# Newton's Method
# Approximates the root of a function "f" with first derivative
# "fp" given some inital guess "p0" within a tolerance of "tol" 
# for a maximum of "N" iterations.
function newtons(f, fp, p0, tol, N)
    p = p0

    # Prints header of output table and initial guess results
    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
    "Error")
    println("--------------------------------------------\
     ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, f(p0))


    for i in 1:N
        func = f(p0)
        deriv = fp(p0)
        
        if deriv == 0 
            # Method fails if root guess is a critical
            return @error("f'(p0) = 0", p0)
        else
            # Newton's method iteration to find next root guess
            p = p0 - f(p0)/deriv
        end
        
        # Estimate error using difference between successive root guesses
        error = p == 0 ? abs(p-p0) : abs(p-p0)/abs(p)
        # Print current results in output table
        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, func, error)

        # Output root guess once estimated error falls within tol
        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end
        
        # Update value of previous root guess
        p0 = p
    end

    return @error("method failed after N interations", N)
end

# Modified Newton's Method
# Approximates the root of a function "f" with first derivative
# "fp" and second derivative "fpp" given some inital guess "p0" 
# within a tolerance of "tol" for a maximum of "N" iterations.
function newtons(f, fp, fpp, p0, tol, N)
    p = p0

    # Prints header of output table and initial guess results
    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
    "Error")
    println("--------------------------------------------\
     ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, f(p0))


    for i in 1:N
        func = f(p0)
        deriv = fp(p0)
        deriv2 = fpp(p0)
        
        # Modified Newton's method uses second derivative values
        # to accelerate convergence
        p = p0 - (func*deriv)/(deriv^2 - func*deriv2)
        
        error = p == 0 ? abs(p-p0) : abs(p-p0)/abs(p)
        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, func, error)

        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end
        
        p0 = p
    end

    return @error("method failed after N interations", N)
end


# Secant Method
# Approximates the root of a function "f" given some inital 
# guesses "p0" and "p1" within a tolerance of "tol" for a 
# maximum of "N" iterations.
#
# "falsepos" is a boolean variable that if set to true will 
# force the algorithm to impement the false position variation 
# of the secant method.
function secantmethod(f, p0, p1, tol, N, falsepos = false)
    q0 = f(p0)
    q1 = f(p1)

    if falsepos && sign(q0)*sign(q1) > 0
        return @error("initial guesses q0 = f(p0) and \
                        q1 = f(p1) must be opposite in sign if \
                        false position is to be used.", q0, q1)
    end

    # Prints header of output table and initial guess results
    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
            "Error")
    println("--------------------------------------------\
             ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, q0)
    @printf("%d\t%.8e\t  %.8e\t%.8e\n", 1, p1, q1, abs(p1-p0))

    for i in 2:N
        # Use secant method to update root guess
        p = p1 - (p1 - p0)*q1/(q1 - q0)
        q = f(p)

        error = p == 0 ? abs(p-p1) : abs(p-p1)/abs(p)

        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, q, error)

        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end

        # If using false position method, p0/q0 guess values will
        # only be updated as needed to ensure that q0 and q1 have opposite
        # sign
        if !falsepos || sign(q)*sign(q1) < 0
            p0 = p1
            q0 = q1
        end 
        p1 = p
        q1 = q
    end

    @error("method failed after N interations", N)
end


# Steffensen's Method
# Approximates the fixed point of some function "g" given
# an inital guess "p0" within a tolerance of "tol" for a 
# maximum of "N" iterations.
function steffensen(g, p0, tol, N)
    # Prints header of output table and initial guess results
    println("\nn", "\t", "p_n", "\t\t", "  g(p_n)", "\t\t", 
            "Error")
    println("--------------------------------------------\
             ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, g(p0))

    for i in 1:N
        p1 = g(p0)
        p2 = g(p1)

        # Use Steffensen's method to update root guess value
        p = p0 - ((p1 - p0)^2) / (p2 - 2*p1 + p0)

        error = p == 0 ? abs(p-p0) : abs(p-p0)/abs(p)

        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, g(p), error)

        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end

        p0 = p

    end

    @error("method failed after N interations", N)
end


# Horner's Method
# Used to quickly evaluate a polynomial with coefficients "a" and
# its first derivative at the point "x0"
#
# Note: "a" must be ordered from coefficients of HIGHEST degree
# to coefficients of LOWEST degree:
# a = [a_n, a_{n-1}, a_{n-2}, ..., a_1, a_0]
#
# Returns a tuple in the form of (P(x0), P'(x0)), where P is the
# polynomial whose coefficients are given by "a"
function horners(a, x0)
    n = length(a) - 1

    bk = a[1] 
    ck = a[1]

    for i in 2:n
        bk = x0*bk + a[i]
        ck = x0*ck + bk
    end

    bk = x0*bk + a[n+1]

    return bk, ck
end
