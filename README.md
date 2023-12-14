# lexmin

Calculation of lexicographically smallest isomorphic model for a given one.

# building

On a Linux machine, the `configure` script should download minisat and cadical
and prepare the `build` folder.

Run `./configure -h`  to see options for building  (for instance a switch to cadical).

In a nutshell, it should be enough to do:
```
     ./configure && cd build && make
```

# usage

Currently the program only support algebras comprising a single binary function, e.g. semigroup.
By default a set of algebras in the GAP format are expetected on the input.
This means that input is a list of list of lists (each algebra is a list of
rows). See the `examples` folder for example inputs.

Run `./mlex -h` for more options.
