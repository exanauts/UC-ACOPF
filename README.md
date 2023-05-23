# UC-ACOPF
A collection of UC-ACOPF models

To initiate a CPU run of the UCMP model, do

```julia runner.jl case9 6```

The first argument is the name of case. Currently, we have data for `case9` and `case118`. The second argument is the number of time periods. Recommended values are 6, 12, 24. 

To initiate a GPU run, do 

```julia runner.jl case9 6 use_gpu```
