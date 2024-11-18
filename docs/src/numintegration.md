# Numerical Integration

```@contents
```

## Composite Methods

```@docs
simpson(f::Function, a, b, n=2)
trapezoid(f::Function, a, b, n=1)
gaussianquad(f, a, b, n; sub=1, nwlist=nothing, checklist=true)
```

## Adaptive Methods

```@docs
adaptivequad(f, a, b, tol, N, base=true)
adaptivegaussquad(f, a, b, n, tol, submax, base=true)
```

```@index
```