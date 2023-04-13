#########################################################
# Interpolate
# File of functions used for creating interpolating
# polynomials for data sets.

# Author: Masen Pitts
# Last Modified: 11/05/2022 (MM/DD/YYYY)
# Version: 1.1
#########################################################

# Neville's Method
# Evaluates the Lagrange interpolating polynomial  at "x_eval" that agrees with some
# function at points "x" with corresponding values of "f"
function neville(x, f, x_eval)
    len = 0

    # Ensure that input arrays are of the same size
    if length(x) != length(f)
        return println("Error: Input vectors x and f must be the same size")
    else
        len = length(f)
    end

    # Initialize first column of Neville table with f(x) values
    Q = zeros((len, len))
    Q[:,1] = f

    # Format output table
    printer = "i\tx_i"

    for k in 1:len
        printer = string(printer, "\t\tQ_i$(k-1)")
    end

    println(printer)

    println("--------------------------------------------\
     --------------------------------------")

    @printf("0\t%.7f\t%.7f\n", x[1], Q[1,1])

    # Display iterated rows of Neville table
    for i in 2:len
        printer = @sprintf("%d\t%.7f", i-1, x[i])
        for j in 2:i
            Q[i,j] = ( (x_eval - x[i-j+1])*Q[i,j-1] - (x_eval - x[i])*Q[i-1,j-1] ) / (x[i] - x[i-j+1])
            printer = string(printer, @sprintf("\t%.7f", Q[i,j-1]))
        end

        printer = string(printer, @sprintf("\t%.7f", Q[i,i]))
        println(printer)
    end

    @printf("\nP(%.2f) = %.7f\n", x_eval, Q[len,len])

    # Approximate value is last entry in Neville table
    return Q[len,len]
end

# Newton's Divided Differences
# Evaluates the Lagrange interpolating polynomial  at "x_eval" that agrees with some
# function at points "x" with corresponding values of "f" using the method
# of Newton's Divided Differences.
function newtondivdiff(x, f, x_eval=0)
    len = 0

    # Ensure that input arrays are of the same size
    if length(x) != length(f)
        return println("Error: Input vectors x and f must be the same size")
    else
        len = length(f)
    end

    # Initialize first column of divided difference table with f(x) values
    F = zeros((len,len))
    F[:,1] = f

    # Format output table
    printer = "i\tx_i"

    for k in 1:len
        printer = string(printer, "\t\tF_i$(k-1)")
    end

    println(printer)

    println("--------------------------------------------\
     --------------------------------------")

    @printf("0\t%.8f\t%.8f\n", x[1], F[1,1])

    # Display iterated rows of divided difference table
    for i in 2:len
        printer = @sprintf("%d\t%.8f", i-1, x[i])
        for j in 2:i
            F[i,j] = (F[i,j-1] - F[i-1,j-1]) / (x[i] - x[i-j+1])
            printer = string(printer, @sprintf("\t%.8f", F[i,j-1]))
        end

        printer = string(printer, @sprintf("\t%.8f", F[i,i]))
        println(printer)
    end

    # Approximate f(x_eval) using Newton's Divided Difference Formula
    P = F[1,1]

    for i in 2:len
        prod = 1

        for j in 1:i-1
            prod *= x_eval - x[j]
        end
        
        P += F[i,i]*prod
    end

    @printf("\nP(%.2f) = %.8f\n", x_eval, P)

    return P
end

# Hermite Interpolation
# Evaluates the Hermite interpolating polynomial  at "x_eval" that agrees with some
# function at points "x" with corresponding values of "f" first derivatives "fp"
# using the method of Newton's Divided Differences.
function hermite(x, f, fp, x_eval)
    len = 0

    # Ensure that input arrays are of the same size
    if length(x) != length(f)
        return println("Error: Input vectors x and f must be the same size")
    elseif length(f) != length(fp) 
        return println("Error: Input vectors f and fp must be the same size")
    elseif length(x) != length(fp)
        return println("Error: Input vectors x and fp must be the same size")
    else
        len = length(x)
    end

    # Format output table
    printer = "i\tz_i"

    for k in 1:2*len
        printer = string(printer, "\t\tQ_i$(k-1)")
    end

    println(printer)

    println("--------------------------------------------\
        --------------------------------------")

    # Determine coefficients of Hermite interpolating polynomial given by the
    # diagonal elements of Q
    Q = zeros((2*len, 2*len))
    z = zeros(2*len)

    for i in 1:len
        z[2*i-1] = x[i]
        z[2*i] = x[i]
        Q[2*i-1,1] = f[i]
        Q[2*i,1] = f[i]
        Q[2*i,2] = fp[i]

        if i != 1
            Q[2*i-1,2] = (Q[2*i-1,1] - Q[2*i-2,1]) / (z[2*i-1] - z[2*i-2])
        end
    end

    for i in 3:2*len
        for j in 3:i
            Q[i,j] = (Q[i,j-1] - Q[i-1,j-1]) / (z[i] - z[i-j+1])
        end
    end

    @printf("0\t%.7f\t%.7f\n", z[1], Q[1,1])

    for i in 2:2*len
        printer = @sprintf("%d\t%.7f", i-1, z[i])

        for j in 1:i
            printer = string(printer, @sprintf("\t%.7f", Q[i,j]))
        end

        println(printer)
    end


    # Approximate value of f(x_eval) using Hermite interpolating polynomial
    H = Q[1,1]

    for i in 2:2*len
        prod = 1

        for j in 1:i-1
            prod *= x_eval - z[j]
        end

        H += Q[i,i]*prod
    end

    @printf("\nH(%.2f) = %.7f\n", x_eval, H)

    return H
end

# Natural Cubic Spline
# Constructs a natural cubic spline interpolant S for the function with
# given values "f" at "x"
function natcubicspline(x, f)
    len = 0

    # Ensure that input arrays are of the same size
    if length(x) != length(f)
        return println("Error: Input vectors x and f must be the same size")
    else
        len = length(f)
    end

    # Solve the linear system for the cubic spline coefficients at each
    # sub-interval
    ###############
    a = f
    h = zeros(len-1)
    alpha = zeros(len-1)

    for i in 1:len-1
        h[i] = x[i+1] - x[i]
    end

    for i in 2:len-1
        alpha[i] = (3/h[i])*(a[i+1] - a[i]) - (3/h[i-1])*(a[i] - a[i-1])
    end

    l = zeros(len)
    mu = zeros(len-1)
    z = zeros(len)

    l[1] = 1
    mu[1] = 0
    z[1] = 0

    for i in 2:len-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i]
    end

    c = zeros(len)
    b = zeros(len-1)
    d = zeros(len-1)
    l[len] = 1
    z[len] = 0
    c[len] = 0

    for j in len-1:-1:1
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    ###############

    # Print cubic spline function
    println("S(x) =")
    for i in 1:len-1
        @printf("[%d] %.5f + %.5f(x - %.1f) + %.5f(x - %.1f)^2 \
                    + %.5f(x - %.1f)^3,\n\tfor x in [%.1f, %.1f]\n", 
                    i, a[i], b[i], x[i], c[i], x[i], d[i], x[i], x[i], x[i+1])
    end

    # List of anonymous functions comprised of the cubic splines at each sub-interval
    S = [s -> a[i] + b[i]*(s - x[i]) + c[i]*(s - x[i])^2 + 
                d[i]*(s - x[i])^3 for i in 1:len-1]

    # Return anonymous function that evaluates the input argument 
    # on the proper interval
    S_x = s -> 
        begin
            i_eval = 1

            for i in 1:len-1
                if s >= x[i] && s <= x[i+1]
                    i_eval = i
                end
            end

            S[i_eval](s)
        end

    return S_x
end

# Clamped Cubic Spline
# Constructs a natural cubic spline interpolant S for the function with
# given values "f" at "x" and endpoint first derivatives "fp0" and "fpn"
function clampcubicspline(x, f, fp0, fpn)
    len = 0

    # Ensure that input arrays are of the same size
    if length(x) != length(f)
        return println("Error: Input vectors x and f must be the same size")
    else
        len = length(f)
    end

    # Solve the linear system for the cubic spline coefficients at each
    # sub-interval
    ###############
    a = f
    h = [x[i+1] - x[i] for i in 1:len-1]
    alpha = zeros(len)

    alpha[1] = 3*(a[2] - a[1])/h[1] - 3*fp0
    alpha[len] = 3*fpn - 3*(a[len] - a[len-1])/h[len-1]

    for i in 2:len-1
        alpha[i] = (3/h[i])*(a[i+1] - a[i]) - (3/h[i-1])*(a[i] - a[i-1])
    end

    l = zeros(len)
    mu = zeros(len-1)
    z = zeros(len)

    l[1] = 2*h[1]
    mu[1] = 0.5
    z[1] = alpha[1]/l[1]

    for i in 2:len-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i]
    end

    c = zeros(len)
    b = zeros(len-1)
    d = zeros(len-1)
    l[len] = h[len-1]*(2 - mu[len-1])
    z[len] = (alpha[len] - h[len-1]*z[len-1])/l[len]
    c[len] = z[len]

    for j in len-1:-1:1
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    ###############

    # Print cubic spline function
    println("S(x) =")
    for i in 1:len-1
        @printf("[%d] %.5f + %.5f(x - %.1f) + %.5f(x - %.1f)^2 \
                    + %.5f(x - %.1f)^3,\n\tfor x in [%.1f, %.1f]\n", 
                    i, a[i], b[i], x[i], c[i], x[i], d[i], x[i], x[i], x[i+1])
    end

    # List of anonymous functions comprised of the cubic splines at 
    # each sub-interval
    S = [s -> a[i] + b[i]*(s - x[i]) + c[i]*(s - x[i])^2 + 
                d[i]*(s - x[i])^3 for i in 1:len-1]

    # Return anonymous function that evaluates the input argument
    # on the proper interval
    S_x = s -> 
        begin
            i_eval = 1

            for i in 1:len-1
                if s >= x[i] && s <= x[i+1]
                    i_eval = i
                end
            end

            S[i_eval](s)
        end

    return S_x
end

# Very cool example
# using Plots

# duckx = [.9,1.3,1.9,2.1,2.6,3.0,3.9,4.4,4.7,5.0,6.0,7.0,8.0,9.2,10.5,11.3,11.6,12.,12.6,13.,13.3];
# duckf = [1.3,1.5,1.85,2.1,2.6,2.7,2.4,2.15,2.05,2.1,2.25,2.3,2.25,1.95,1.4,0.9,0.7,0.6,0.5,0.4,0.25];
# duck = natcubicspline(duckx, duckf);
# xGrid = 0.9:0.1:13.3; yVals = [duck(x) for x in xGrid];
# plot(xGrid, yVals, xlims=(0,14), ylims=(-6,4))
# plot!(duckx, duckf, seriestype = :scatter)
