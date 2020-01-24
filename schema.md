


Want to have a unified table, largely cnosisting exactly of what pops out of
Andy's code.

NOTE: All s, mu, nu, etc. values are in the analytic normalization.



Main table
----------

| Field | type | description |
|----------|    ------     | ----- |
label | string | ???
origin | string | an LMFDB url from which we have derived this data
conductor | numeric
is_motivic | bool
motivic_weight | smallint
degree | smallint
mus | numeric[][] | the shifts of gamma_r (double array for Maass forms, etc.). NOTE: must make decision regarding mu = [0,1] vs nu = [1]
nus | numeric[][] | the shifts of gamma_c (double array for Maass forms, etc.)
primitive | bool |
root_number_arg | numeric[] | in some range (as in arb)
primitive_central_character | int[] | pair of integers in Conrey labelling scheme
self_dual | bool |
conjugate_label | string | probably the Lhash
Lhash | string | hash of lowest zero of L-function
leading_term | numeric[] | first nonzero coefficient of the Taylor expansion centered at critical point
values | numeric[] | list of pairs (s_0, L(s_0)), where L(s_0) is stored as arb


## removed
 - symmetry_type


Other table
-----------
is_algebraic | bool | whether the coeffs and euler_factors will be stored as ints
Euler_factors |
bad_euler_factors |
coeffs |


There is a question to be resolved about whether or not we store the coefficients of L-functions that
are algebraic numbers in some number field as algebraic numbers, or as complexes that "forget" where
they come from. In one sense, L-functions are first class citizens and should be considered on their
own right. One could conceive of an exercise where we compute all L-functions with a particular sort
of transformation law up to some conductor bound; in this case, it might be of interest whether these
are algebraic or not.

Closely related is whether we should store Galois orbits of L-functions, in a structure roughly
analogous to how modular forms are stored.

Per Drew, it would be nice to store and show algebraic integers when possible. This should be done
for generic Galois-orbit L-function; and particular L-function embeddings should know their embedding.
Further, each element of the orbit should know its orbit and link to it.

