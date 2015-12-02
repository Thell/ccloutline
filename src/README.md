---
title: "Notes"
---

## Flow

- R NumericMatrix convertable by RcppArmadillo to a C++ Armadillo matrix.
  - 8-way connected border states
    - connected component labeling
    - lattice (grid) graph edge vertice pairs
      - Strongly connected component eulerian cycles vertice trails
        - Cartesian coordinate eulerian cycle paths
          - R List of numeric matrices indicating label, component and path.

## Notes

~~~

Consider that the input matrix is `n_row` rows by `n_col` columns,
each unique `r x c` pair is a cell and each cell has four border
edges between four corner points. These points define the lattice.
That lattice (as a matrix) is indexed such that the first index
represents the vertex index of the Northwest corner point of the first
cell in the input matrix.
The index values in the lattice proceed down the column until the
`n_row + 1`th row. The last vertex index for the lattice matrix
represents the Southeast corner point of the last cell from the input.

There are a total of `(n_row+1)(n_col+1)` vertices in the lattice.


The strongly connected component cycle extraction routine which
uses the input edges, a lookup table with space for each incidental
vertex, two position tracking stacks and the resulting cycles
container.

Armadillo uses an [`arma::uword`][1] for indexing.

- `uword` is a typedef for an `unsigned integer type`; it is used
  for matrix indices as well as all internal counters and loops.
- when using the new C++11 / C++14 standards, the default width
  is 64 bits on 64-bit platforms (as of Armadillo 5.000)


Unfortunately, it is defined as `unsigned long long` even when the
c++ fundamental `size_t` type is defined as `unsigned long` creating
a conversion mismatch. On top of this, Rcpp does not provide a
builtin wrap for `size_t` or `arma::uword`, so conversion to `double`
is needed before passing the results back to R.


The largest addressable 64bit `size_t` addressable value is
$2^{64} - 1$, therefore the maximum number of input values must be
limited to $2^{64} - 1 - n\_row -n\_col$. Since the input _must_ be
at least two dimensions greater than 1



[1]:http://arma.sourceforge.net/docs.html#uword
