To run the experiments, do

```bash
julia --project=. run_experiments.jl > output.log 2>&1
```



1) rotations
2) check intersecting sets and non intersecting sets --> 
    a) generate polytope A randomly, this in high enough dimensions will probably look like a ball
    b) generate another polytope B (maybe like a pyramid or something pointy) and take a vertex x from inside A, possibly close to the "surface"
        of the ball-looking polytope: the pointy vertex from B and x must coincide    
3) print lmo characteristics on logs
4) iterations, dual gap, function over time, time
5) ask mathieu about BPCG differences w(block coordinate)
6) comparison with BC: you can find a point in the intersection using our method: see COLT paper