#########################################################
# Iterative Linear Algebra
# File of functions using iterative numerical methods
# to solve linear systems.

# Author: Masen Pitts
# Last Modified: 12/06/2022 (MM/DD/YYYY)
# Version: 1.0
#
# Primary contents:
#   jacobi      23
#   gaussseidel 100
#########################################################

# Include various linear algebra helper methods
include("lin_algebra.jl")

# Jacobi Iterative Method
# Solves the matrix equation
#   A*x = b
# within tolerance "tol" given initial guess vector "x0"
# and maximum number of iterations "Niter," where "A" is 
# an N x N matrix.
function jacobi(A, b, x0, tol, Niter)
    N = length(b)

    A = typecheckfloat(A)

    # Check that sizes of vectors and matrices are appropriate
    if N != size(A)[1]
        return @error("Input matrix A and input vector b must have compatible shapes.")
    elseif N != length(x0)
        return @error("Input vector b and input vector x0 must be of the same size.")
    elseif size(A)[1] != size(A)[2]
        return @error("Input matrix A must be square (dimensions N x N).")
    end

    x0 = x0[:] # Make local variable "x0" so that input vector is not mutated
    x = zeros(N) # Vector where solution values are stored

    # Format output table
    printStr = @sprintf("k\t")
    for i in 1:N
        printStr = string(printStr, @sprintf("x_%d\t\t", i))
    end
    printStr = string(printStr, "Error")

    println(printStr)

    # Print initial approximation
    printStr = @sprintf("0\t")
    for i in 1:N
        printStr = string(printStr, @sprintf("%.6f\t", x0[i]))
    end
    printStr = string(printStr, "n/a")

    println(printStr)

    # Fail the method if algorithm runs N times without converging within "tol"
    for k in 1:Niter

        # Use Jacobi Method to find new elements of solution vector
        for i in 1:N
            x[i] = (-innerprod(A[i,:],x0) + A[i,i]*x0[i] + b[i] )/A[i,i]
        end

        # Display results of current iteration
        printStr = @sprintf("%d\t", k)
        for i in 1:N
            printStr = string(printStr, @sprintf("%.6f\t", x[i]))
        end

        # Calculate error of current iteration
        error = norminfvector(x-x0)/norminfvector(x)

        printStr = string(printStr, @sprintf("%.6e", error))

        println(printStr)

        # Check to see if tolerance has been reached
        if error < tol
            print("\n")
            for i in 1:N
                @printf("x_%d = %.8f\n", i, x[i])
            end

            return x
        end

        x0 = x[:] # Set new "x0" to old "x" values
    end

    return @error("Exceeded maximum number of iterations.", Niter)
end

# Gauss-Seidel Iterative Method
# Solves the matrix equation
#   A*x = b
# within tolerance "tol" given initial guess vector "x0"
# and maximum number of iterations "Niter," where "A" is 
# an N x N matrix.
function gaussseidel(A, b, x0, tol, Niter)
    N = length(b)

    A = typecheckfloat(A)

    # Check that sizes of vectors and matrices are appropriate
    if N != size(A)[1]
        return @error("Input matrix A and input vector b must have compatible shapes.")
    elseif N != length(x0)
        return @error("Input vector b and input vector x0 must be of the same size.")
    elseif size(A)[1] != size(A)[2]
        return @error("Input matrix A must be square (dimensions N x N).")
    end

    x0 = x0[:] # Make local variable "x0" so that input vector is not mutated
    x = zeros(N) # Vector where solution values are stored

    # Format output table
    printStr = @sprintf("k\t")
    for i in 1:N
        printStr = string(printStr, @sprintf("x_%d\t\t", i))
    end
    printStr = string(printStr, "Error")

    println(printStr)

    # Print initial approximation
    printStr = @sprintf("0\t")
    for i in 1:N
        printStr = string(printStr, @sprintf("%.6f\t", x0[i]))
    end
    printStr = string(printStr, "n/a")

    println(printStr)

    # Fail the method if algorithm runs N times without converging within "tol"
    for k in 1:Niter

        # Use Gauss-Seidel Method to find new elements of solution vector
        for i in 1:N
            if i == 1
                x[i] = (-innerprod(A[i,i+1:N],x0[i+1:N]) + b[i])/A[i,i]
            else
                x[i] = (-innerprod(A[i,1:i-1],x[1:i-1])-innerprod(A[i,i+1:N],x0[i+1:N]) + b[i])/A[i,i]
            end
        end

        # Display results of current iteration
        printStr = @sprintf("%d\t", k)
        for i in 1:N
            printStr = string(printStr, @sprintf("%.6f\t", x[i]))
        end

        # Calculate error of current iteration
        error = norminfvector(x-x0)/norminfvector(x)

        printStr = string(printStr, @sprintf("%.6e", error))

        println(printStr)

        # Check to see if tolerance has been reached
        if error < tol
            print("\n")
            for i in 1:N
                @printf("x_%d = %.8f\n", i, x[i])
            end

            return x
        end

        x0 = x[:] # Set new "x0" to old "x" values
    end

    return @error("Exceeded maximum number of iterations.", Niter)
end