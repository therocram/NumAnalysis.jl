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
include("direct_lin_algebra.jl")

function initprint(x0, token1, token2)
    # Format output table
    printStr = string("m", token2)
    for i in eachindex(x0)
        printStr = string(printStr, @sprintf("x_%d", i), token1)
    end
    printStr = string(printStr, "muh", token1, "Error")

    println(printStr)

    # Print initial approximation
    printStr = string("0", token2)
    for i in eachindex(x0)
        printStr = string(printStr, @sprintf("%.6f", x0[i]), token2)
    end
    printStr = string(printStr, "n/a", token1,"n/a")

    println(printStr)
end

function printeigen(k, x, muhat, error, token2)
    # Print currnet approximation
    printStr = string(k,token2)
    for i in eachindex(x)
        printStr = string(printStr, @sprintf("%.6f", x[i]), token2)
    end
    
    printStr = string(printStr, @sprintf("%.6f", muhat), token2)

    # Print out error of current iteration
    printStr = string(printStr, @sprintf("%.6e", error))

    println(printStr)
end

# Power Method
# Uses the power method to approximate the dominant eigenvalue
# of square matrix (dimensions: N x N) "A" within tolerance "tol",
# given initial eigenvector guess "x". Will run for a maximum of "Niter"
# iterations before failing. The "usesymm" parameter, which defaults to true,
# tells the function to use the symmetric power method if the input matrix "A"
# is symmetric
#
# Note: This function uses Aitken's delta^2 procedure to accelerate convergence
function eigenpower(A, x, tol, Niter, usesymm=true)
    N = size(A)[1]

    # Essentially converts x to a vector type to allow compatible 
    # shape with matrix A
    x = x[:]

    # Check matrix and vector for compatible sizes
    if N != size(A)[2]
        return @error("Input matrix A must be square", A)
    elseif N != length(x)
        return @error("Input matrix A and input vector x must have \
        compatible sizes", A, x)
    end

    # If matrix is symmetric then default to using Symmetric Power Method
    # unless instructed otherwise
    if usesymm && issymmetric(A)
        return eigenpowersym(A, x, N, tol, Niter)
    end

    # Initialize variables used in Aitken's procedure
    mu0 = 0
    mu1 = 0

    # Scale "x" by its first dominant component
    p = domcomponent(x)
    x /= x[p]

    # Special whitespace strings used for table formatting
    token1 = "        "
    token2 = "   "

    # Format output table and print initial approximation
    initprint(x, token1, token2)
    
    # Main algorithm starts. Runs until tolerance is satisfied or until 
    # Niter iterations have been done
    for k in 1:Niter
        # Iterate next eigenvalue estimate
        y = A*x
        mu = y[p]
        muhat = mu0 - (mu1 - mu0)^2/(mu - 2*mu1 + mu0)

        # Update value of p to be index of dominant component of y
        p = domcomponent(y)

        # Quit if 0 is an eigenvalue
        if y[p] == 0
            return @warn("A has the eigenvalue 0, select a new vector \
            x and restart", x)
        end

        # Calculate error and iterate next eigenvector estimate
        yyp = y/y[p]
        error = norminfvector(x - yyp)
        x = yyp

        # Print out results of current iteration
        printeigen(k, x, muhat, error, token2)

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

    return @error("Exceeded maximum number of iterations.", Niter)

end

# Symmetric Power Method
# Uses the symmetric power method to approximate the dominant eigenvalue
# of symmetric square matrix (dimensions: N x N) "A" within tolerance "tol",
# given initial eigenvector guess "x". Will run for a maximum of "Niter"
# iterations before failing. 
# This function is expected to only be used as a helper method for the 
# "eigenpower" function.
#
# Note: This function uses Aitken's delta^2 procedure to accelerate convergence
function eigenpowersym(A, x, N, tol, Niter)
    # Normalize input vector x
    x /= norm(x)

    # Initialize variables used in Aitken's procedure
    mu0 = 0
    mu1 = 0

    # Special whitespace strings used for table formatting
    token1 = "        "
    token2 = "   "

    # Format output table and print initial approximation
    initprint(x0, token1, token2)

    # Main algorithm starts. Runs until tolerance is satisfied or untile Niter
    # iterations have been done
    for k in 1:Niter
        # Iterate next eigenvalue estimate
        y = A*x
        mu = innerprod(x,y)
        muhat = mu0 - (mu1 - mu0)^2/(mu - 2*mu1 + mu0)

        # Quit if 0 is an eigenvalue
        if norm(y) == 0
            return @warn("A has the eigenvalue 0, select a new vector x and restart", x)
        end

        # Calculate error and iterate next eigenvector estimate
        ynorm = y/norm(y)
        error = norm(x - ynorm)
        x = ynorm

        printeigen(k, x, muhat, error, token2)

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
        mu1 = mu
    end
    
    return @error("Exceeded maximum number of iterations.", Niter)

end

# Inverse Power Method
# Uses the inverse power method to approximate an eigenvalue
# and associated eigenvector of input matrix "A" within tolerance "tol",
# given initial eigenvector guess "x". Will run for a maximum of "Niter"
# iterations before failing.
#
# Note: This function uses Aitken's delta^2 procedure to accelerate convergence
function eigeninvpower(A, x, tol, Niter)
    N = size(A)[1]

    # Essentially converts x to a vector type to allow compatible 
    # shape with matrix A
    x = x[:]

    # Check matrix and vector for compatible sizes
    if N != size(A)[2]
        return @error("Input matrix A must be square", A)
    elseif N != length(x)
        return @error("Input matrix A and input vector x must have compatible sizes", A, x)
    end

    # Initialize variables used in Aitken's procedure
    mu0 = 0
    mu1 = 0

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

    # Calculate eigenvalue guess and scale x by its dominant component
    q = innerprod(x,A*x)/innerprod(x,x)
    p = domcomponent(x)
    x /= x[p]
    
    # Use LU factorization to factor the matrix 
    #   A - qI
    # This speeds up the algorithm since the program
    # must solve the system
    #   (A - qI)y = x
    # at each step
    factoredAqI = matrixfactorLU(A-q*IMatrix(N))

    # If (A - qI) is singular then q must be an eigenvalue of A
    if factoredAqI === nothing
        return @warn("q is an eigenvalue. Select a new vector x and restart", q, x)
    end

    # If (A - qI) is nonsingular L and U are defined as the lower
    # and upper triangular matrices respectively with elements such that
    # LU = (A - qI)
    L, U = factoredAqI

    # Main algorithm starts. Runs until tolerance is satisfied or untile Niter
    # iterations have been done
    for k in 1:Niter
        # Solve the system
        # (A - qI)y = x
        # for y using the factor matrices L and U
        y = systemsolveLU(L, U, x)

        # Next eigenvalue estimate
        mu = y[p]
        muhat = mu0 - (mu1 - mu0)^2/(mu - 2*mu1 + mu0)

        # Next eigenvector esitmate
        p = domcomponent(y)
        yyp = y/y[p]
        error = norminfvector(x - yyp) # Estimate error
        x = yyp

        # Display results of current iteration
        printStr = @sprintf("%d\t", k)
        for i in 1:N
            printStr = string(printStr, @sprintf("%.6f\t", x[i]))
        end
        printStr = string(printStr, @sprintf("%.6f\t", 
                            muhat == 0 ? muhat : 1/muhat + q))

        # Print out error of current iteration
        printStr = string(printStr, @sprintf("%.6e", error))

        println(printStr)  

        # Print results and end process if tolerance is satisfied
        if error < tol && k >= 4
            muhat = 1/muhat + q

            print("\n")
            for i in 1:N
                @printf("x_%d = %.8f\n", i, x[i])
            end
            @printf("muhat = %.8f\n", muhat)

            return muhat, x
        end

        # Update Aitken variables
        mu0 = mu1
        mu1 = mu        
    end

    return @error("Exceeded maximum number of iterations.", Niter)

end

# Housholder's Method
# Returns a symmetric tridiagonal matrix similar to symmetric input matrix
# A using Housholder's Method
function housholder(A)
    # Verify that A is symmetric
    if !issymmetric(A)
        return @error("Input matrix A must be symmetric", A)
    end

    A = typecheckfloat(A)

    N = size(A)[1]

    v = zeros(N)
    u = zeros(N)
    z = zeros(N)

    # Implement algorithm to update elements of A as described in
    # pg. 598 of Numerical Analysis Burden and Faires, 9th Edition,
    # using vectorized operation where possible
    for k in 1:N-2
        q = innerprod(A[k+1:N,k], A[k+1:N,k])

        Ak1k = A[k+1,k]

        if Ak1k == 0
            alpha = -sqrt(q)
        else
            alpha = -sqrt(q)*Ak1k/abs(Ak1k)
        end

        rsq = alpha*alpha - alpha*Ak1k

        v[k] = 0
        v[k+1] = Ak1k - alpha

        v[k+2:N] = A[k+2:N,k]

        for j in k:N
            u[j] = innerprod(A[j,k+1:N], v[k+1:N])/rsq
        end

        prod = innerprod(v[k+1:N], u[k+1:N])

        z[k:N] = u[k:N] - (prod/(2*rsq))*v[k:N]

        k1N1 = k+1:N-1

        for i in k1N1
            i1N = i+1:N

            A[i1N,i] = A[i1N,i] - v[i]*z[i1N] - v[i1N]*z[i]
            A[i,i1N] = A[i1N,i]

            A[i,i] = A[i,i] - 2*v[i]*z[i]
        end

        A[N,N] = A[N,N] - 2*v[N]*z[N]

        A[k,k+2:N] = A[k+2:N,k] .= 0

        A[k+1,k] = A[k+1,k] - v[k+1]*z[k]
        A[k,k+1] = A[k+1,k]
    end

    return A
end