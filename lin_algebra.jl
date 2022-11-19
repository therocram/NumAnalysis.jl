##########################################################
# Linear Algebra
# File of helper functions using direct numerical methods
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

    return transpose(P)*L, U
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


#######################################################################
# Linear Algebra Helper Methods
# A collection of functions performing useful operations
# on vectors (1D arrays) and matrices (2D arrays)
#######################################################################

# Perform gauss jordan elimination on matrix "A"
# NOTE: If checktype = "false" method will permute "A"
function gaussjordanelim(A, checktype = true)
    if checktype
        A = typecheckfloat(A)
    end

    N = size(A)[1]

    # Iterate over columns of matrix and perform
    # elimination
    for i in 1:N-1
        # Search for first nonzero value in the column
        k = 0
        for j in i:N
            if A[j,i] != 0
                k = j
                break
            end
        end

        # Swap rows if desired pivot is in another column
        if k == 0
            return @error("Input matrix \"A\" is singular.")
        elseif k != i
            swaprow(A, i, k) # Swap rows as necessary
        end

        # Eliminate column elements beneath A_ii
        reducecol(A, i)
    end

    return A
end

# Returns true if "A" is a positive definite matrix. Returns false otherwise.
function isposdef(A)
    # Must be symmetric to be positive definite
    if !issymmetric(A)
        return false
    end

    A = typecheckfloat(A)

    N = size(A)[1]

    # Iterate over columns of matrix and perform
    # elimination
    for i in 1:N-1
        # "A" is positive definite if an only if Gaussian elimination can be
        # without row interchanges can be performed with all pivot elements positive
        if A[i,i] <= 0
            return false
        end
        reducecol(A, i)
    end

    return true
end

# Returns true if "A" is a symmetric matrix. Returns false otherwise.
function issymmetric(A)
    if A == transpose(A)
        return true
    else
        return false
    end
end

# Returns N x N Identity matrix
function IMatrix(N)
    I = zeros(N,N)

    for i in 1:N
        for j in 1:N
            I[i,j] = ==(i,j)
        end
    end

    return I
end

# Checks data type of matrix and converts its entries to floating point
# numbers as necessary 
function typecheckfloat(A)
    if isa(A, Matrix{Int})
        A = float.(A)
        @warn("Converting integer elements to floating point elements.\n\
               Avoid inputing matrix of pure integers if possible.")
    end

    return A[:,:]
end

# Swaps rows "n" and "m" of the matrix "A"
function swaprow(A, n, m)
    for i in 1:size(A)[2]
        A[n,i], A[m,i] = A[m,i], A[n,i]
    end
end

# Swaps columns "n" and "m" of the matrix "A"
function swapcol(A, n, m)
    for i in 1:size(A)[1]
        A[i,n], A[i,m] = A[i,m], A[i,n]
    end
end

# Performs the row operation:
#   (n + w*m) -> n
# on the matrix "A"
function addrow(A, n, m, w)
    for i in 1:size(A)[2]
        A[n,i] = A[n,i] + w*A[m,i]
    end
end

# Eliminates the column elements beneath A_ii
# NOTE: Assumes that column has been properly pivoted/scaled
# as desired
function reducecol(A, i)
    for j in i+1:size(A)[1]
        addrow(A, j, i, -A[j,i]/A[i,i])
    end
end

# Create augmented N x N+1 matrix "Ab" using N x N matrix "A"
# and N x 1 (or 1 x N) vector "b"
function augmentmatrix(A, b)
    N = size(A)[1]

    if N != size(A)[2]
        return @error("Input matrix \"A\" must be square (dimensions N x N).")
    end

    Ab = zeros((N,N+1))
    
    # Fill all entries of "Ab" except last column with entries of "A"
    # Then fill last column with enetries of "b"
    Ab[:,1:N] = A
    Ab[:,N+1] = b

    return Ab
end

# Performs forward subsitution on reduced N x N+1 augmented
# matrix and returns list of solution values "x" in the order
# x = [x1, x2, ..., xN]
function forwardsub(A, printx=false)
    N = size(A)[1]

    if (A[N,N] == 0) 
        return @error("No unique solution exists", A) 
    end

    # Vector storing solution values
    x = zeros(N)

    # Perform forwards substitution on reduced matrix
    x[1] = A[1,N+1]/A[1,1]
    for i in 2:N
        x[i] = (A[i,N+1] - innerprod(A[i,1:i-1], x[1:i-1])) / A[i,i]
    end

    # Format a printed result if desired
    if printx
        for j in 1:N
            @printf("x%d = %.7f\n", j, x[j])
        end
    end

    return x
end

# Performs backward subsitution on reduced N x N+1 augmented
# matrix and returns list of solution values "x" in the order
# x = [x1, x2, ..., xN]
function backsub(A, printx=false)
    N = size(A)[1]

    if (A[N,N] == 0) 
        return @error("No unique solution exists", A) 
    end

    # Vector storing solution values
    x = zeros(N)

    # Perform backwards substitution on reduced matrix
    x[N] = A[N,N+1]/A[N,N]
    for i in N-1:-1:1
        x[i] = (A[i,N+1] - innerprod(A[i,i+1:N], x[i+1:N])) / A[i,i]
    end

    # Format a printed result if desired
    if printx
        # line1 = @sprintf("i\t")
        # line2 = @sprintf("x_i\t")

        # for j in 1:N
        #     line1 = string(line1, @sprintf("%d\t\t", j))
        #     line2 = string(line2, @sprintf("%.4f\t\t", x[j]))
        # end

        # println(line1,"\n",line2)

        for j in 1:N
            @printf("x%d = %.7f\n", j, x[j])
        end
    end

    return x
end

# Evaluates the inner product (row_vector*col_vector product)
# of "a" and "b" over real scalars
function innerprod(a, b)
    if length(a) != length(b)
        return @error("\"a\" and \"b\" must be vectors of the same length")
    end

    sum = 0
    for i in eachindex(a)
        sum += a[i]*b[i]
    end 

    return sum
end

# Returns the maximum element in a 1D list
function maxlist(v)
    max = v[1]

    for i in 2:lastindex(v)
        if v[i] > max
            max = v[i]
        end
    end

    return max
end
