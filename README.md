## ccloutline

Connected component labeling, euler cycle extraction and outline
point paths.

----

Processing of a 10,000x10,000 matrix with a random binary population
takes __400mb of input__ data and generates __7.3gb of output__ in
~28s on a i7-4700MQ with 24gb memory single threaded.

The connected component labeling took ~3s for 656,972 individual
components. The eulerian path extraction (200,012,368 edges) represented
an additional 9 seconds after generating component edges. The largest
consumer of time is the final list of matrices structure wrapping.

Total memory consumed while processing peaked at just under 19gb.

Rather than a linked list or adjacency matrix for edges, each edge
in the lattice structure is directed and each vertex is stored with
its outgoing neighbors in a lookup array of tuples; taking advantage
of the fact that there is at most 2 degrees for any vertex.

See:

- CCLREMSP - Algorithm 4 _A New Parallel Algorithm for Two-Pass Connected
Component Labeling_, Gupta et al.

- Cheriyan-Mehlhorn-Gabow SCC Algorithm - An Extended Experimental Evaluation of SCC (Gabow's vs Kosaraju's) based on Adjacency List
 
---------------------

### Example

While the examples just show binary input, any values within a numeric
matrix that are comparable for equality will work.

input:

~~~
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0
0 1 1 1 1 1 1 1 1 0 0 1 1 1 1 0 0
0 0 0 1 1 1 1 0 0 0 1 1 1 1 0 0 0
0 0 1 1 1 1 0 0 0 1 1 1 0 0 1 1 0
0 1 1 1 0 0 1 1 0 0 0 1 1 1 0 0 0
0 0 1 1 0 0 0 0 0 1 1 0 0 0 1 1 0
0 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
~~~

labels and component populations:

~~~
> ccl_labels(sample_mat)
[[1]]
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
 [1,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1     1     1     1
 [2,]    1    1    2    2    1    1    2    2    1     1     3     3     1     1     3     3     1
 [3,]    1    2    2    2    2    2    2    2    2     1     1     3     3     3     3     1     1
 [4,]    1    1    1    2    2    2    2    1    1     1     3     3     3     3     1     1     1
 [5,]    1    1    2    2    2    2    1    1    1     3     3     3     1     1     3     3     1
 [6,]    1    2    2    2    1    1    2    2    1     1     1     3     3     3     1     1     1
 [7,]    1    1    2    2    1    1    1    1    1     3     3     1     1     1     3     3     1
 [8,]    1    1    1    1    1    1    3    3    3     3     1     1     3     3     3     3     1
 [9,]    1    1    1    1    1    1    1    1    1     1     1     1     1     1     1     1     1

[[2]]
[1] 94 27 32
~~~

input:

~~~
0 1 0 1 0
1 0 1 0 1
0 1 0 1 0
1 0 1 0 1
0 1 0 1 0
~~~

Eulerian cycles _of the cell vertices_ in column major order:

~~~
> ccl_cycles(sample_mat)
> ccl_cycles(tst128x128cb2[1:5,1:5])
[[1]]
[[1]][[1]]
 [1]  0  6  7 13 12 18 19 13 14 20 19 25 24 30 31 25 26 20 21 27 26 32 33 27 28 34 35
[28] 29 28 22 23 17 16 10 11  5  4 10  9  3  2  8  9 15 16 22 21 15 14  8  7  1  0


[[2]]
[[2]][[1]]
 [1]  1  7  6 12 13  7  8 14 13 19 18 24 25 19 20 14 15 21 20 26 25 31 32 26 27 21 22
[28] 28 27 33 34 28 29 23 22 16 17 11 10  4  3  9 10 16 15  9  8  2  1
~~~

outlines:

input:

~~~
0 0 0
0 1 0
0 0 0
~~~

output:

~~~
[[1]]
[[1]][[1]]
      label ccl x y
 [1,]     0   1 0 0
 [2,]     0   1 0 1
 [3,]     0   1 0 2
 [4,]     0   1 0 3
 [5,]     0   1 1 3
 [6,]     0   1 2 3
 [7,]     0   1 3 3
 [8,]     0   1 3 2
 [9,]     0   1 3 1
[10,]     0   1 3 0
[11,]     0   1 2 0
[12,]     0   1 1 0
[13,]     0   1 0 0

[[1]][[2]]
     label ccl x y
[1,]     0   2 1 1
[2,]     0   2 2 1
[3,]     0   2 2 2
[4,]     0   2 1 2
[5,]     0   2 1 1


[[2]]
[[2]][[1]]
     label ccl x y
[1,]     1   1 1 1
[2,]     1   1 1 2
[3,]     1   1 2 2
[4,]     1   1 2 1
[5,]     1   1 1 1
~~~

The outlines can be fed into a plotting function to create borders
for each component.
