# Calculus
# File of helper functions numerically approximating derivatives
# and integrals

# Composite Simpson's Rule.
# Approximates the integral Int(f(x), a < x < b) using
# Simpson's Rule on n subintervals.
function simpson(f, a, b, n)
    h = (b - a)/n

    y0 = f(a) + f(b)
    y1 = 0; y2 = 0

    for i in 1:n-1
        y = a + i*h

        if i % 2 == 0
            y2 += f(y)
        else
            y1 += f(y)
        end
    end

    return h*(y0 + 2*y2 + 4*y1)/3
end