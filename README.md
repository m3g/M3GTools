## Distribution function

Given a list of values, for example, 
```
v = rand(1000)
```
creates a distribution function of the values, with:
```
using M3GTools
x, y = distribution(v)
```
Parameters `vmin` and `vmax` define the limits of the distribution.
Parameter `step` defines the width (in units of `v`) of the histogram
count. By default 100 bins will be created and `step=(vmax-vmin)/100`.
If `step=5`, the step will be 5 times greater.

In general, one would use
```
using M3GTools
x, y = distribution(v,step=5)
```
and adjust the `step` parameter to obtain a nice distribution function, considering that
at each `x` the `y` value will contain the probability of finding a value of `v`
within `[x-step,x+step]`, where `step=5*(vmax-vmin)/100` (which only sums up to 1.0 if
`step = (vmax-vmin)/nbins`, that is, if there is no overlap between the counts. 


