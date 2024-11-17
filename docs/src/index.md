# Your worst nightmare

```@docs
simpson(f::Function, a, b, n=2)
trapezoid(f::Function, a, b, n=1)
adaptivequad(f, a, b, tol, N, base=true)
gaussianquad(f, a, b, n; sub=1, nwlist=nothing, checklist=true)
adaptivegaussquad(f, a, b, n, tol, submax, base=true)
ivpsolve(f, a, b, y0, N; solver=rungekuttaO4, y=nothing, printres=false)
```