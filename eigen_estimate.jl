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
# given initial eigenvector guess "x" and initial eigenvalue guess "q."
# If "q" is not explicitly provided then an inital guess will be generated inside
# the method.
# Will run for a maximum of "Niter" iterations before failing.
#
# Note: This function uses Aitken's delta^2 procedure to accelerate convergence
function eigeninvpower(A, x, tol, Niter, q = nothing)
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

    # Special whitespace strings used for table formatting
    token1 = "        "
    token2 = "   "

    # Format output table
    initprint(x, token1, token2)

    # Come up with initial eigenvalue guess q if one is not explicitly provided
    if q === nothing
        q = innerprod(x,A*x)/innerprod(x,x)
    end
    # Calculate eigenvalue guess and scale x by its dominant component
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
        printeigen(k, x, muhat, error, token2)

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

    # Ensure that A has a compatible data type
    A = typecheckfloat(A)

    N = size(A)[1]

    # Initialize some helper arrays
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
            j = i+1:N

            A[j,i] = A[j,i] - v[i]*z[j] - v[j]*z[i]
            A[i,j] = A[j,i]

            A[i,i] = A[i,i] - 2*v[i]*z[i]
        end

        A[N,N] = A[N,N] - 2*v[N]*z[N]

        A[k,k+2:N] = A[k+2:N,k] .= 0

        A[k+1,k] = A[k+1,k] - v[k+1]*z[k]
        A[k,k+1] = A[k+1,k]
    end

    return A
end

# Returns similar upper hessenberg matrix of arbitrary input matrix A
function hessenberg(A)
    # Ensure that A has a compatible data type
    A = typecheckfloat(A)

    N = size(A)[1]

    # Initialize some helper arrays
    v = zeros(N)
    u = zeros(N)
    z = zeros(N)
    y = zeros(N)

    # Implement algorithm to update elements of A as described in
    # pg. 598-600 of Numerical Analysis Burden and Faires, 9th Edition,
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

        for j in 1:N
            u[j] = innerprod(A[j,k+1:N], v[k+1:N])/rsq
            y[j] = innerprod(A[k+1:N,j], v[k+1:N])/rsq
        end

        prod = innerprod(v[k+1:N], u[k+1:N])

        z = u - (prod/rsq)*v

        k1N = k+1:N

        for i in k+1:N
            j = 1:k

            A[j,i] = A[j,i] - z[j]*v[i]
            A[i,j] = A[i,j] - y[j]v[i]

            j = k1N

            A[j,i] = A[j,i] - z[j]*v[i] - y[i]*v[j]
        end
    end

    return A
end

# Formatting for printing eigenvalue approximations and appending
# eigenvalues to ongoing list in QR algorithm
function seteigen(eigenvals, eigen, M)
    append!(eigenvals, eigen)
    @printf("\u03bb_%d = %.8f\n", M - length(eigenvals) + 1, eigen)
end

# QR Algorithm
# Obtains the eigenvalues of the symmetric, tridiagonal input matrix "A"
# using the QR algorithim within the given tolerance "tol" for a maximum of 
# "Niter" iterations
#
# Note: Shifting, which accelerates the algorithm's convergence, is implemented
# by default but can be disabled by setting "shifting" to false
function qreigen(A, tol, Niter, shifting=true)
    # Verify that A is symmetric and tridiagonal
    if !issymmetric(A)
        return @error("Input matrix A must be symmetric", A)
    elseif !istridiag(A)
        return @error("Input matrix A must be tridiagonal", A)
    end

    N = size(A)[1]
    M = N # Constant holding the original size of the input matrix 

    # Populate arrays storing diagonal and off-diagaonal elements of A
    a = zeros(N)
    b = zeros(N)

    # Helper arrays used in main algorithm
    f = zeros(N)
    x = zeros(N)
    y = zeros(N)
    z = zeros(N)
    g = zeros(N)
    s = zeros(N)
    q = zeros(N)
    r = zeros(N)

    # Used to store and return eigenvalue approximations
    eigenvals = zeros(0)
    #m = 0 # keeps track of the number of eigenvalues generated

    for i in 1:N
        a[i] = A[i,i]
        if i != N  b[i+1] = A[i,i+1] end
    end

    #display(buildtridiag(a, b, b))

    shift = 0 # This shift factor whill remain zero unless shifting = true
    #simga = 0

    # Run main algorithm for a maximum of N iterations
    for k in 1:Niter
        # Checks for eigenvalues within tolerance in last row of current matrix
        if abs(b[N]) <= tol 
            seteigen(eigenvals, a[N] + shift, M)
            N -= 1 # Decrease size of current matrix by 1 in each dimension
        end

        # Checks for eigenvalues within tolerance in second row of current matrix
        if abs(b[2]) <= tol
            seteigen(eigenvals, a[1] + shift, M)
            N -= 1 # Decrease size of current matrix by 1 in each dimension
            a[1:N] = a[2:N+1] # Update values of a and b to reflect change in
            b[2:N] = b[3:N+1] # matrix size
        end

        # Ends algorithm and returns eigenvalues if matrix has
        # been reduced to nothing
        if N == 0
            return eigenvals
        end

        # Ends algorithm and returns eigenvalues if matrix has
        # been reduced to one element
        if N == 1
            seteigen(eigenvals, a[1] + shift, M)
            return eigenvals
        end

        # Splits up A into two seperate matrices and returns eigenvalues
        # if off diagonal elements within tolerances are detected that are not 
        # at the beginning or end of the matrix
        for j in 3:N-1
            if b[j] <= tol
                println("Split matrix into")
                display(buildtridiag(a[1:j-1], b[2:j-1], b[2:j-1]))
                println("and")
                display(buildtridiag(a[j:N], b[j+1:N], b[j+1:N]))

                return eigenvals
            end
        end

        # Compute shifting factor if shifting = true
        if shifting
            e = -(a[N-1] + a[N])
            c = a[N]*a[N-1] - b[N]*b[N]
            d = sqrt(e*e - 4*c)

            if e > 0
                mu1 = -2*c/(e + d)
                mu2 = -(e + d)/2
            else
                mu1 = (d - e)/2
                mu2 = 2*c/(d - e)
            end
            #println(mu1, " ", mu2)

            # Stop algorithm and return eigenvalues if matrix has been reduced to
            # 2 x 2 matrix
            if N == 2
                seteigen(eigenvals, mu1 + shift, M)
                seteigen(eigenvals, mu2 + shift, M)

                return eigenvals
            end

            # Set sigma to closest value to a[N] out of mu1 and mu2
            sigma = abs(mu1 - a[N]) < abs(mu2 - a[N]) ? mu1 : mu2

            # Accumulate shift amount
            shift += sigma

            for j in 1:N
                f[j] = a[j] - sigma
            end
        else
            for j in 1:N
                f[j] = a[j]
            end
        end
        
        ##### Compute upper triangular R matrix
        x[1] = f[1]
        y[1] = b[2]

        for j in 2:N
            z[j-1] = sqrt(x[j-1]*x[j-1] + b[j]*b[j])
            g[j] = x[j-1]/z[j-1]
            s[j] = b[j]/z[j-1]
            q[j-1] = g[j]*y[j-1] + s[j]*f[j]
            x[j] = -s[j]*y[j-1] + g[j]*f[j]

            if j != N
                # r[j-1] = s[j]*b[j+1]
                y[j] = g[j]*b[j+1]
            end
        end
        #####

        ##### Compute next iteration of A
        z[N] = x[N]
        a[1] = s[2]*q[1] + g[2]*z[1]
        b[2] = s[2]*z[2]

        for j in 2:N-1
            a[j] = s[j+1]*q[j] + g[j]*g[j+1]*z[j]
            b[j+1] = s[j+1]*z[j+1]
        end

        a[N] = g[N]*z[N]
        #####
    end

    @printf("A^(%d):\n", Niter + 1)
    display(buildtridiag(a, b[2:M], b[2:M]))
    return @error("Exceeded maximum number of iterations.", Niter)    

end