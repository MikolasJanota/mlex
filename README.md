# lexmin

The tool calculates the lexicographically smallest isomorphic model for a given
model (algebra). The model is given as a Cayley (multiplication) table. The calculation
is done using a SAT solver as described in [1].

# Building

On Linux , the `configure` script should download minisat and cadical
and prepare the `build` folder.

Run `./configure -h`  to see options for building  (for instance a switch to cadical).

In a nutshell, it should be enough to do:
```
     ./configure && cd build && make
```

# Usage

Currently the program only support algebras comprising a single binary function, e.g. a semigroup.
By default, a set of algebras in the GAP format are expected on the input.
This means that input is a list of list of lists (each algebra is a list of
rows). The MACE4 format is also supported. See the `examples` folder for example inputs.

Run `./mlex -h` for more options.

# References

[1] SAT-based Techniques for Lexicographically Smallest Finite Models in AAAI 2024.
M. Janota, C. Chow, J. Araújo, M. Codish, P. Vojtěchovský
