# lexmin

Calculation of lexicographically smallest isomorphic model.

# building

`setup-release.sh` script should download minisat and cadical and build in the `build` folder on a linux machine

To switch to cadical, do `cd build && cmake -DUSE_CADICAL=ON && make` 

# usage

By default gap format is expetected on the input. Issu `./mlex -h` for more options.
