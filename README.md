# Digital Convolution
A common operation in scientific computing is digital convolution, by which each element in a multi-dimensional grid is replaced by a weighted sum of its neighbors. This has applications in graphics, in which such operations are done to both blur and sharpen images, and in numerical simulations, in which this may be a single step in a Poisson solver.

# OpenMPI
OpenMPI is a popular implementation of a message passing interface, here used on TACC to increase performance of our digital convolution algorithm.

# Results
Scaling was overall close to linear. There is an initial cost to distributing data across tasks and synchronizing computational loads, but comparing N=2 to N=32 sees a linear increase in performance. See attached plots

### Example output (Ran with 8 nodes): 
Num elements in an array col                   : 16386
Num inner elements in an array col             : 16384
Num elements in an array                       : 268500996
Num inner elements in an array                 : 268435456
Num bytes used by an array (GB)                : 1.074004
Num elements below threshold in array y        : 2990
Fraction of elements below threshold in array y: 0.000011
Num elements below threshold in array x        : 26849775
Fraction of elements below threshold in array x: 0.099999
Threshold value:                               : 0.100000

Timing summary (Seconds).
WALL: Alloc-X :  0.000054
WALL: Alloc-Y :  0.000390
WALL: Init-X  :  18.261909
WALL: Smooth-Y:  0.409883
WALL: Count-Y :  0.021454
WALL: Count-X :  0.020173

