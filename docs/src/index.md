# Your worst nightmare

```@docs
gaussianquad(f, a, b, n; sub=1, nwlist=nothing, checklist=true)
adaptivequad(f, a, b, tol, N, base=true)
adaptivegaussquad(f, a, b, n, tol, submax, base=true)
ivpsolve(f, a, b, y0, N; solver=rungekuttaO4, y=nothing, printres=false)
```