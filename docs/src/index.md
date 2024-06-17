# Your worst nightmare

```@docs
gaussianquad(f, a, b, n; sub=1)
adaptivequad(f, a, b, tol, N, base=true)
ivpsolve(f, a, b, y0, N; solver=rungekuttaO4, y=nothing, printres=false)
```