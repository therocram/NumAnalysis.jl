########################################################
# Eigen Estimate
# File of functions using various numerical methods
# to solve for the eigenvalues and eigenvectors of linear
# systems.

# Author: Masen Pitts
# Last Modified: 1/19/2023 (MM/DD/YYYY)
# Version: 1.0
##########################################################

using Printf

# Include various linear algebra helper methods
include("lin_algebra.jl")

# Power Method
# Uses the power method to approximate the dominant eigenvalue
# of square matrix (dimensions: N x N) "A" within tolerance "tol",
# given initial eigenvector guess "x". Will run for a maximum of "Niter"
# iterations before failing
#
# Note: This function uses Aitken's delta^2 procedure to accelerate convergence
function eigenpower(A, x, tol, Niter)
    N = size(A)[1]

    # Check matrix and vector for compatible sizes
    if N != size(A)[2]
        return @error("Input matrix A must be square", A)
    elseif N != length(x)
        return @error("Input matrix A and input vector x must have compatible sizes", A, x)
    end

    # Initialize variables used in Aitken's procedure
    mu0 = 0
    mu1 = 0

    # Scale "x" by its first dominant component
    p = domcomponent(x)
    x /= x[p]

    # Format output table
    printStr = @sprintf("m\t")
    for i in 1:N
        printStr = string(printStr, @sprintf("x_%d\t\t", i))
    end
    printStr = string(printStr, "muhat\t\tError")

    println(printStr)

    # Print initial approximation
    printStr = @sprintf("0\t")
    for i in 1:N
        printStr = string(printStr, @sprintf("%.6f\t", x[i]))
    end
    printStr = string(printStr, "n/a\t\tn/a")

    println(printStr)
    
    # Fail the method if algorithm runs N times without converging within "tol"
    for k in 1:Niter
        # Iterate next eigenvalue estimate
        y = A*x
        mu = y[p]
        muhat = mu0 - (mu1 - mu0)^2/(mu - 2*mu1 + mu0)

        # Display results of current iteration
        printStr = @sprintf("%d\t", k)
        for i in 1:N
            printStr = string(printStr, @sprintf("%.6f\t", x[i]))
        end
        printStr = string(printStr, @sprintf("%.6f\t", muhat))

        # Update value of p to be index of dominant component of y
        p = domcomponent(y)

        # Quit if 0 is an eigenvalue
        if y[p] == 0
            return @warn("A has the eigenvalue 0, select a new vector x and restart", x)
        end

        # Calculate error and iterate next eigenvector estimate
        yyp = y/y[p]
        error = norminfvector(x - yyp)
        x = yyp

        # Print out error of current iteration
        printStr = string(printStr, @sprintf("%.6e", error))

        println(printStr)

        # Check to see if tolerance has been reached
        if error < tol && k >= 4
            print("\n")
            for i in 1:N
                @printf("x_%d = %.8f\n", i, x[i])
            end
            @printf("muhat = %.8f\n", muhat)

            return muhat, x
        end

        # Update Aitken's procedure variables
        mu0 = mu1
        mu1 =  mu
    end
end

function eigenpowersym()
end