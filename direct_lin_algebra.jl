########################################################
# Direct Linear Algebra
# File of functions using direct numerical methods
# to solve linear systems.

# Author: Masen Pitts
# Last Modified: 11/21/2022 (MM/DD/YYYY)
# Version: 1.0
#
# Primary contents:
#   gaussbacksub        25
#   gausspivotpartial   48
#   gausspivotscaled    93
#   matrixfactorLU      148
#   matrixfactorPtLU    207
#   matrixfactorLDLt    267
#   matrixfactorLLt     304
#   systemsolveLU       344
#########################################################

using Printf

# Include various linear algebra helper methods
include("lin_algebra.jl")

# Gaussian Elimination with Backwards Substitution
# Directly solves a linear system by performing
# Gaussian elimination on a provided N x N+1 augmented
# matrix "A" describing the system.
#
# Returns solutions "x" in the order
# x = [x1, x2, ..., xN]
function gaussbacksub(A)
    A = typecheckfloat(A)

    N = size(A)[1] # Number of rows in A

    if size(A)[2] != N+1
        return @error("Augmented matrix must have shape n x n+1", A)
    end

    # Iterate over columns of augmented matrix and perform
    # elimination
    gaussjordanelim(A, false)

    backsub(A, true)
end

# Gaussian Elimination with Partial Pivoting
# Directly solves a linear system by performing
# Gaussian elimination with partial pivoting on a provided 
# N x N+1 augmented matrix "A" describing the system.
#
# Returns solutions "x" in the order
# x = [x1, x2, ..., xN]
function gausspivotpartial(A)
    A = typecheckfloat(A)

    N = size(A)[1] # Number of rows in A

    if size(A)[2] != N+1
        return @error("Augmented matrix must have shape n x n+1", A)
    end

    # Iterate over columns of augmented matrix and perform
    # elimination
    for i in 1:N-1
        # Search column pivot
        k = 0
        maxA = 0
        for j in i:N
            aji = A[j,i]
            if abs(aji) > maxA
                k = j
                maxA = aji
            end
        end

        # Swap rows if desired pivot is in another column
        if k == 0
            return @error("No unique solution exists", A) 
        elseif k != i
            swaprow(A, i, k) # Swap rows as necessary
        end

        # Eliminate column elements beneath A_ii
        reducecol(A, i)

    end

    return backsub(A, true)
end

# Gaussian Elimination with Scaled Partial Pivoting
# Directly solves a linear system by performing
# Gaussian elimination with scaled partial pivoting on a provided 
# N x N+1 augmented matrix "A" describing the system.
#
# Returns solutions "x" in the order
# x = [x1, x2, ..., xN]
function gausspivotscaled(A)
    A = typecheckfloat(A)

    N = size(A)[1] # Number of rows in A

    if size(A)[2] != N+1
        return @error("Augmented matrix must have shape n x n+1", A)
    end

    # Determine scaling factors for each row
    s = zeros(N)
    for i in 1:N
        s[i] = maxlist(abs.(A[i,1:N])) # Scaling factor is abs value of maximum entry in row

        if s[i] == 0
            @error("No unique solution exists", A) # Only happens if row is all zeros
        end
    end

    # Iterate over columns of augmented matrix and perform
    # elimination
    for i in 1:N-1
        # Search for column pivot with largest scaled value
        k = 0
        maxA = 0
        for j in i:N
            aji = abs(A[j,i])/s[j]
            if aji > maxA
                k = j
                maxA = aji
            end
        end

        # Swap rows if desired pivot is in another column
        if k == 0
            return @error("No unique solution exists", A) 
        elseif k != i
            swaprow(A, i, k)
            s[i], s[k] = s[k], s[i]
        end

        # Eliminate column elements beneath A_ii
        reducecol(A, i)
    end

    return backsub(A, true)
end

# LU Factorization
# Factors the N x N matrix "A" into the product of a lower
# triangular matrix "L" with 1 for all of its diagonal entries
# and an upper triangular matrix "U."
# A = L*U
# NOTE: Uses Doolittle's Method of LU factorization
function matrixfactorLU(A)
    N = size(A)[1] # Number of columns and rows in "A"

    # "A" must be a square matrix
    if N != size(A)[2]
        return @error("Input matrix A must be square matrix (dimensions N x N).")
    end

    A = typecheckfloat(A)

    L = zeros((N,N)) # Lower triangular matrix
    U = zeros((N,N)) # Upper triangular matrix

    # Use Doolittle's Method of LU factorization
    L[1,1] = 1; U[1,1] = A[1,1]

    if L[1,1]*U[1,1] == 0
        return @error("Factorization impossible. Input matrix has a zero diagonal element.")
    end

    # Set first row of "U" and first column of "L"
    for j in 2:N
        U[1,j] = A[1,j]/L[1,1]
        L[j,1] = A[j,1]/U[1,1]
    end

    for i in 2:N-1
        # Calculate diagonal elements of "L" and "U" using Doolittle's Method
        L[i,i] = 1; U[i,i] = A[i,i] - innerprod(L[i,1:i-1], U[1:i-1,i])

        if L[i,i]*U[i,i] == 0
            return @error("Factorization impossible. \
                           Upper triangular matrix has a zero diagonal element", U, L)
        end

        # Set ith row of "U" and ith column of "L"
        for j in i+1:N
            U[i,j] = (A[i,j] - innerprod(L[i,1:i-1], U[1:i-1,j]))/L[i,i]
            L[j,i] = (A[j,i] - innerprod(L[j,1:i-1], U[1:i-1,i]))/U[i,i]
        end
    end

    # Set last diagonal elements of "L" and "U" using Doolittle's Method
    L[N,N] = 1; U[N,N] = A[N,N] - innerprod(L[N,1:N-1], U[1:N-1,N])

    if L[N,N]*U[N,N] == 0
        @warn("Input Matrix \"A\" is singular.")
    end

    return L, U
end

# PtLU Factorization
# Factors the N x N matrix "A" into the product of the transpose of a permutation matrix "P" 
# with a lower triangular matrix "L" with 1 for all of its diagonal entries
# and an upper triangular matrix "U."
# A = (Pt*L)*U
# NOTE: Uses Doolittle's Method of LU factorization
function matrixfactorPtLU(A)
    N = size(A)[1]

    # "A" must be a square matrix
    if N != size(A)[2]
        return @error("Input matrix A must be square matrix (dimensions N x N).")
    end

    A = typecheckfloat(A)
    U = A[:,:] # Upper triangular matrix
    L = IMatrix(N) # Lower triangular matrix (Doolittle method)
    P = IMatrix(N) # Permutation matrix

    # Iterate over columns of matrix and perform
    # elimination to find set elements of "U" and "P"
    for i in 1:N-1
        # Search for first nonzero value in the column
        k = 0
        for j in i:N
            if U[j,i] != 0
                k = j
                break
            end
        end

        # Swap rows if desired pivot is in another column
        if k == 0
            return @error("factorization is impossible.", U)
        elseif k != i
            swaprow(U, i, k) # Swap rows as necessary
            swaprow(P, i, k) # Any swap done with U must also be done with P
        end

        # Eliminate column elements beneath U_ii
        reducecol(U, i)
    end

    PA = P*A

    # Perform Gaussian elimination on PA and determine the elements
    # of "L"
    for i in 1:N-1
        for j in i+1:N
            w = PA[j,i]/PA[i,i]
            L[j,i] = w
            addrow(PA, j, i, -w)
        end
    end

    return P, L, U
end

# LDLt Factorization
# Factors the N x N positive definite matrix "A" into the product of a lower triangular matrix "L" 
# with 1 for all of its diagonal entries with a diagonal matrix "D" and "Lt," the transpose of "L."
# A = L*D*Lt
function matrixfactorLDLt(A)
    # Check if A is positive definite
    if !isposdef(A)
        return @error("Input matrix \"A\" must be positive definite.")
    end

    A = typecheckfloat(A)

    N = size(A)[1]

    L = IMatrix(N) # Lower triangular matrix with ones on diagonal
    D = zeros(N,N) # Diagonal matrix
    v = zeros(N-1) # Vector used in factorization

    # Fill the entries of "L" and "D"
    for i in 1:N
        for j in 1:i-1
            v[j] = L[i,j]*D[j,j]
        end

        D[i,i] = A[i,i] - innerprod(L[i,1:i-1],v[1:i-1])

        for j in i+1:N
            L[j,i] = (A[j,i] - innerprod(L[j,1:i-1],v[1:i-1]))/D[i,i]
        end
    end

    L, D, transpose(L)
end

# LLt Factorization
# Factors the N x N positive definite matrix "A" into the product of a lower triangular matrix "L" 
# and its transpose "Lt."
# A = L*Lt
# NOTE: Uses Cholesky's Method
function matrixfactorLLt(A)
    # Check if A is positive definite
    if !isposdef(A)
        return @error("Input matrix \"A\" must be positive definite.")
    end

    A = typecheckfloat(A)

    N = size(A)[1]
    L = zeros((N,N)) # Lower triangular matrix

    L[1,1] = sqrt(A[1,1]) # Set first element of "L" (Cholesky's Method)

    # Set first column of "L"
    for j in 2:N
        L[j,1] = A[j,1]/L[1,1]
    end

    for i in 2:N-1
        # Set diagonal elements of "L"
        L[i,i] = sqrt(A[i,i] - innerprod(L[i,1:i-1],L[i,1:i-1]))

        # Set corresponding column of "L"
        for j in i+1:N
            L[j,i] = (A[j,i] - innerprod(L[j,1:i-1],L[i,1:i-1]))/L[i,i]
        end
    end

    # Set last diagonal entry of "L"
    L[N,N] = sqrt(A[N,N] - innerprod(L[N,1:N-1],L[N,1:N-1]))

    return L, transpose(L)
end

# System Solve LU
# Solves the system 
# Ax = b
# for "x," where A = L*U has been factored into the product of a lower triangular
# matrix "L" and an upper triangular matrix "U."
function systemsolveLU(L, U, b)
    N = size(L)[1]

    if N != size(L)[2]
        return @error("Input matrices must be square (dimensions N x N).", L)
    elseif N != size(U)[1] || N != size(U)[2]
        return @error("Input matrices \"L\" and \"U\" must be of the same size.", L)
    end

    L = typecheckfloat(L)
    U = typecheckfloat(U)

    # Solve the equation L*y = b using forward subsitution
    Ly = augmentmatrix(L, b)
    y = forwardsub(Ly)
    

    # Solve the equation U*x = y using backwards substitution
    Ux = augmentmatrix(U, y)
    x = backsub(Ux, true)

    return x
end
