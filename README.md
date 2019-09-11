# Distribution function

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

