##########################################################
# Root Finder
# File of helper functions used for finding the roots and
# fixed points of real functions of one variable

# Author: Masen Pitts
# Last Modified: 11/05/2022 (MM/DD/YYYY)
# Version: 1.1
#########################################################


# Bisection Method
# Implements the bisection method on the interval [a,b]
# to approximate the root of a function "f" within some tolerance
# "tol" for a maximum of "N" iterations.
# "printi" is a boolean variable that determines whether
# the function will print the total number of iterations
# when it returns a value.
function bisecmethod(f, a, b, tol, N)
    if f(a) == 0
        return a
    elseif f(b) == 0
        return b
    end

    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
    "Error")
    println("--------------------------------------------\
     ---------------")

    if f(a)f(b) > 0
        println("Error: Function values at endpoints must \
                have different signs")
        return nothing
    end

    lastp = a

    for i = 1:N
        p = a + (b-a)/2
        val = f(p)
        error = p == 0 ? abs(p-lastp) : abs(p-lastp)/abs(p) 

        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, val, error)

        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end

        if sign(f(a))*sign(val) > 0
            a = p
        else
            b = p
        end

        lastp = p
    end

    println("\nMethod failed after N interations. N = ", N)
end

# Fixed Point
# Approximates the fixed point of some function "g" given
# an inital guess "p0" within a tolerance of "tol" for a 
# maximum of "N" iterations.
function fixedpoint(g, p0, tol, N)
    for i in 1:N
        p = g(p0)
        if p == 0 ? abs(p-p0) < tol : abs(p-p0)/abs(p) < tol
            if printi println("Required Iterations: ", i) end
            return p
        end

        p0 = p
    end

    println("Method failed after N interations. N = ", N)
end

# Newton's Method
# Approximates the root of a function "f" with first derivative
# "fp" given some inital guess "p0" within a tolerance of "tol" 
# for a maximum of "N" iterations.

using Printf

function newtons(f, fp, p0, tol, N)
    p = p0

    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
    "Error")
    println("--------------------------------------------\
     ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, f(p0))


    for i in 1:N
        func = f(p0)
        deriv = fp(p0)
        
        if deriv == 0 
            return @printf("\nError: f'(%f) = 0", p0)
        else
            p = p0 - f(p0)/deriv
        end
        
        error = p == 0 ? abs(p-p0) : abs(p-p0)/abs(p)
        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, func, error)

        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end
        
        p0 = p
    end

    println("\nMethod failed after N interations. N = ", N)
end

# Modified Newton's Method
# Approximates the root of a function "f" with first derivative
# "fp" and second derivative "fpp" given some inital guess "p0" 
# within a tolerance of "tol" for a maximum of "N" iterations.
function newtons(f, fp, fpp, p0, tol, N)
    p = p0

    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
    "Error")
    println("--------------------------------------------\
     ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, f(p0))


    for i in 1:N
        func = f(p0)
        deriv = fp(p0)
        deriv2 = fpp(p0)
        
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

    println("\nMethod failed after N interations. N = ", N)
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
        return println("Error: Initial guesses f(p0) and \
                        f(p1) must be opposite in sign if \
                        false position is to be used.")
    end

    println("\nn", "\t", "p_n", "\t\t", "  f(p_n)", "\t\t", 
            "Error")
    println("--------------------------------------------\
             ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, q0)
    @printf("%d\t%.8e\t  %.8e\t%.8e\n", 1, p1, q1, abs(p1-p0))

    for i in 2:N
        p = p1 - (p1 - p0)*q1/(q1 - q0)
        q = f(p)

        error = p == 0 ? abs(p-p1) : abs(p-p1)/abs(p)

        @printf("%d\t%.8e\t  %.8e\t%.8e\n", i, p, q, error)

        if error < tol
            println("\nRequired Iterations: ", i)
            println("p = ", p)
            return p
        end

        if !falsepos || sign(q)*sign(q1) < 0
            p0 = p1
            q0 = q1
        end 
        p1 = p
        q1 = q
    end

    println("\nMethod failed after N interations. N = ", N)
end


# Steffensen's Method
# Approximates the fixed point of some function "g" given
# an inital guess "p0" within a tolerance of "tol" for a 
# maximum of "N" iterations.
function steffensen(g, p0, tol, N)
    println("\nn", "\t", "p_n", "\t\t", "  g(p_n)", "\t\t", 
            "Error")
    println("--------------------------------------------\
             ---------------")
    @printf("%d\t%.8e\t  %.8e\tn/a\n", 0, p0, g(p0))

    for i in 1:N
        p1 = g(p0)
        p2 = g(p1)

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

    println("\nMethod failed after N interations. N = ", N)
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
