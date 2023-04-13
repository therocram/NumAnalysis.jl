#########################################################
# Nonlinear Solver
# File of functions used for finding the solutions to
# systems of nonlinear equations
#
# Author: Masen Pitts
# Last Modified: 2/08/2022 (MM/DD/YYYY)
# Version: 1.0
#########################################################

# Prints a formatted output table and initial approximation values "x0"
# for iterative numerical methods using whitespace specifiers "token1"
# and "token2"
function initprint(x0, token1, token2)
    # Format output table
    printStr = string("k", token2)
    for i in eachindex(x0)
        printStr = string(printStr, @sprintf("x_%d", i), token1)
    end
    printStr = string(printStr, "Error")

    println(printStr)

    # Print initial approximation
    printStr = string("0", token2)
    for i in eachindex(x0)
        printStr = string(printStr, @sprintf("%.6f", x0[i]), token2)
    end
    printStr = string(printStr, "n/a")

    println(printStr)
end

# Prints the formatted results "x" of the "k"th iteration of an iterative
# numerical algorithm with error "error" using the whitespace specifier "token2"
function printiter(k, x, error, token2)
    # Print currnet approximation
    printStr = string(k,token2)
    for i in eachindex(x)
        printStr = string(printStr, @sprintf("%.6f", x[i]), token2)
    end

    # Print out error of current iteration
    printStr = string(printStr, @sprintf("%.6e", error))

    println(printStr)
end

# A formatted printing of each element of a list of 
# floating point numbers "x"
function printresults(x)
    print("\n")
    for i in eachindex(x)
        @printf("x_%d = %.8f\n", i, x[i])
    end
end

# Passes the input vector "x" to the transformation "F", which
# is a vector whose entries are functions mapping "x" onto the
# real number line, such that the output is also a vector.
function passto(x, F::Vector{Function}, forcefloat=false)
    if forcefloat
        val = Float64[]
    else
        val = []        
    end

    for i in eachindex(F)
        append!(val, F[i](x))
    end

    return val
end

# Passes the input vector "x" to the matrix transformation "F", which
# is a matrix whose entries are functions mapping "x" onto the
# real number line, such that the output is also a matrix.
function passto(x, F::Matrix{Function}, forcefloat=false)
    if forcefloat
        val = Matrix{Float64}(undef, size(F)...)
    else
        val = Matrix{}(undef, size(F)...)
    end

    for i in eachindex(F)
        val[i] = F[i](x)
    end

    return val
end

# Non-Linear Fixed Point Algorithm
# Solves for the fixed point of input transformation "F" 
# (vector of functions) with initial guess "p0" until error
# within "tol" is reached or until "Niter" iterations have occured.
#
# The boolean parameter "usegaussSeidel", which is true by default,
# enables use of the Gauss-Seidel method to accelerate convergence.
function fixedpointNL(F, p0, tol, Niter, usegaussSeidel=true)
    # Ensure that p0 and F are of compatible sizes
    if length(p0) != length(F)
        @error("Input guess vector p0 must be same length \
        as input transformation F", p0, F)
    end

    # Special whitespace strings used for table formatting
    token1 = "        "
    token2 = "   "

    # Format output table and print initial approximation
    initprint(p0, token1, token2)

    for k in 1:Niter
        # By default uses Gauss-Seidel version of fixed point algorithm
        # unless otherwise specified
        if usegaussSeidel
            p = p0[:]

            for j in eachindex(F)
                p[j] = F[j](p)
            end
        else
            p = passto(p0,F,true)
        end

        # Estimate relative error of current approximation
        # error = norminfvector(p0) == 0 ? norminfvector(p - p0) :
        #         norminfvector(p - p0)/norminfvector(p0)
        error = norminfvector(p - p0)

        # Print results of current iteration
        printiter(k, p, error, token2)

        # Output current approximation if error is within tolerance
        if error < tol
            printresults(p)
            return p
        end

        # Store current p value for use in next approximation
        p0 = p
    end

    @error("Exceeded maximum number of iterations", Niter)
end

# Newton's Method for Nonlinear systems
# Approximates the solution of the nonlinear system
#   F(x) = 0
# Where "F" is the input transformation (vector of functions). 
# Other parameters include the initial guess vector "x", 
# the Jacobian transformation matrix "J" (matrix of functions),
# the desired tolerance "tol", and the maximum number of iterations "Niter"
function newtonsNL(F::Vector{Function}, x, J::Matrix{Function}, tol, Niter)
    # Ensure that x0 and F are of compatible sizes
    if length(x) != length(F)
        @error("Input guess vector x must be same length \
        as input transformation F", x, F)
    elseif length(F) != size(J)[1]
        @error("Input transformation F and input Jacobian J \
        must have compatible shapes", x, F)
    end

    # Special whitespace strings used for table formatting
    token1 = "        "
    token2 = "   "

    # Format output table and print initial approximation
    initprint(x, token1, token2)

    for k in 1:Niter
        # Compute entries of J(x)
        Jx = passto(x,J,true)

        # Compute entries of F(x)
        Fx = passto(x,F,true)

        # Forms the augmented matrix [J(x):-F(x)]
        JFx = augmentmatrix(Jx,-Fx)

        # Solve the linear system J(x)*y = -F(x)
        y = gausspivotscaled(JFx)
        
        # Update solution estimate and calculate error
        x += y
        error = norminfvector(y)

        # Output current iteration to data table
        printiter(k, x, error, token2)

        # Output current estimation if error is within tolerance
        if error < tol
            printresults(x)
            return x
        end
    end

    @error("Exceeded maximum number of iterations", Niter)
end
