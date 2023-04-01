#######################################################################
# Linear Algebra Helper Methods
# A collection of functions performing useful operations
# on vectors (1D arrays) and matrices (2D arrays)

# Author: Masen Pitts
# Last Modified: 2/05/2023 (MM/DD/YYYY)
# Version: 1.4
#######################################################################

# using Printf

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

# Returns true if "A" is a tridiagonal matrix. Returns false otherwise
function istridiag(A)
    # Cutoff value. Many algorithms for determining tridiagonal matrices
    # return very small values in empty matrix entries instead of exactly 0.
    # This cutoff value (on order of 1e-10 to 1e-20) takes this into account by
    # considering any matrix entry below the cutoff value to be zeros
    cutoff = 1e-10

    for i in 1:size(A)[1]
        for j in i+2:size(A)[2]
            if A[i,j] >= cutoff
                return false
            end
        end

        for k in i+2:size(A)[2]
            if A[k,i] >= cutoff
                return false
            end
        end
    end

    return true
end

# Creates a tridiagonal matrix with diagonal elements given
# by input vector "a", subdiagonal elements given by "b",
# and superdiagonal elemetns given by "c"
function buildtridiag(a, b, c)
    N = length(a)

    A = zeros(N,N)

    for i in 1:N
        A[i,i] = a[i]
        if i < N
            A[i+1,i] = b[i]
            A[i,i+1] = c[i]
        end
    end

    return A
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

# Create augmented M x N+1 matrix "Ab" using M x N matrix "A"
# and M x 1 (or 1 x M) vector "b"
function augmentmatrix(A, b)
    sA = size(A)

    # if N != size(A)[2]
    #     return @error("Input matrix \"A\" must be square (dimensions N x N).", A)
    # elseif N != length(b)
    #     return @error("Input matrix \"A\" and input vector \"b\" must have compatible shapes", A, b)
    # end

    if sA[1] != length(b)
        return @error("Input matrix \"A\" and input vector \"b\" must have compatible shapes", A, b)
    end  

    Ab = typeof(A)(undef,sA[1],sA[2]+1)

    # Fill all entries of "Ab" except last column with entries of "A"
    # Then fill last column with enetries of "b"
    Ab[:,1:sA[2]] = A
    Ab[:,sA[2]+1] = b

    return Ab
end

# Performs forward subsitution on reduced N x N+1 augmented
# matrix and returns list of solution values "x" in the order
# x = [x1, x2, ..., xN]
function forwardsub(A, printx=false)
    N = size(A)[1]

    if (A[1,1] == 0) 
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

# Returns the norm (square root of sum of squares) of a vector "v"
function norm(v)
    return sqrt(innerprod(v,v))
end

# Returns the infinity norm of some vector "v"
function norminfvector(v)
    return maxlist(abs.(v))
end

# Returns the index of the first component v_p in "v" such that
#   abs(v_p) = infnorm(v)
function domcomponent(v)
    infnorm = norminfvector(v)

    for i in eachindex(v)
        if abs(v[i]) == infnorm
            return i
        end
    end
end

# Returns the natural infinity norm of some N x N matrix "A"
function norminfmatrix(A)
    N = size(A)[1]

    if N != size(A)[2]
        @error("Input matrix A must be square (dimensions N x N)")
    end 

    sums = [sum(abs.(A[i,:])) for i in 1:N]

    return maxlist(sums)    
end
