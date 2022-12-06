# aSlope

A small library to perform simple slope stability calculations
Code adapted to methods given in NTNU course BYGX2001
- this is by no means a finished projects so check all results carefully


## Instructions
The main.py script contains examples of calculations for a simple geometry.

## Examples
This code can be used to calculate factor of safety for sircular cylindric failure surfaces.

It can handle undrained strength parameters (Su) - here shown for an arbitrary failure surface.
![](https://raw.githubusercontent.com/siggimar/aSlope/main/examples/Single_FS_Su.PNG)


And it can handle drained strength parameters as well.
![](https://raw.githubusercontent.com/siggimar/aSlope/main/examples/Single_FS_a-phi.PNG)

It is possible to identify the most critical failure surface and calculate its factor of safety using a simple grid search.
![](https://raw.githubusercontent.com/siggimar/aSlope/main/examples/ex_2_gridsearch.png)
