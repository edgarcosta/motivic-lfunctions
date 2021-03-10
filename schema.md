
Search table: lfunc_search
--------------------------

| Field | type | description |
|----------|    ------     | ----- |
origin                        |text                | url for the object that was use to generated this data
primitive                     |boolean             | true if L-func is [primitive], we use the second moment in many instances to decide this
conductor                     |numeric             | the [conductor] of L-func
central_character             |text                | the label of the primitive [central character](https://beta.lmfdb.org/knowledge/show/lfunction.central_character)
self_dual                     |boolean             | true if L-func is self-dual (coeff field is totally real)
motivic_weight                |smallint            | the [motivic weight] of the L-func 
Lhash                         |text                | **to remove**
degree                        |smallint            | the [degree] of the L-func
order_of_vanishing            |smallint            | the [analytic rank], the order of vanishing at its central point
algebraic                     |boolean             | if the L-func is [arithmetic], i.e. normalized Dirichlet coefficients are algebraic numbers, conjecturally this is the same as being algebraic
z1                            |numeric             | the [lowest zero]
gamma_factors                 |jsonb               | **to remove**
trace_hash                    |bigint              | linear combination of the a_p between 2^12 and 2^13 reduced mod 2^61-1 as defined in BSSVY, only for rational L-functions
root_angle                    |double precision    | stored between -.5 to .5
prelabel                      |text                | the label without the index
analytic_conductor            |double precision    | the analytic conductor [analytic conductor]
mu_real                       |smallint[]          | the real part (in [0, 1]) of mus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
mu_imag                       |numeric[]           | the imaginary part of mus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
nu_real_doubled               |smallint[]          | the real part of mus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
nu_imag                       |numeric[]           | the imaginary part of mus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
bad_primes                    |bigint[]            | primes dividing the [conductor]
label                         |text                | the [label] of the lfunc, e.g.: 3-1-1.1-r0e3-p4.23p33.33m37.56-0
index                         |smallint            | the last component of the [label]
conductor_radical             |integer             | product of the bad_primes, i.e., the primes dividing the [conductor]
dirichlet_coefficients        |numeric[]           | the [Dirichlet coefficients] for rational L-functions in arithmetic normalisation starting with a_1
euler_factors                 |numeric[]           | **to remove??**
rational                      |boolean             | if the L-func is [algebraic]
euler{p}                      |numeric[]           | arrays of length degree + 1 representing the [euler factors] for p = 2, 3, 5,..., 97
root_analytic_conductor       |double precision    | the [root analytic conductor]
instance_types                |text[]              | representing the keys of the multimap url(type) -> url(instance)
instance_urls                 |text[]              | representing the values of the multimap url(type) -> url(instance)
spectral_label                |text                | the [spectral label], e.g. r0e3-p4.23p33.33m37.56
is_instance_{type}            |boolean             | type in {Artin,BMF,CMF,DIR,ECNF,G2Q,HMF,MaassGL3,MaassGL4,MaassGSp4,NF}

[label]: https://beta.lmfdb.org/knowledge/show/lfunction.label
[spectral label]: https://beta.lmfdb.org/knowledge/show/lfunction.spectral_label
[analytic rank]: https://beta.lmfdb.org/knowledge/show/lfunction.analytic_rank
[lowest zero]: https://beta.lmfdb.org/knowledge/show/lfunction.zeros
[degree]: https://beta.lmfdb.org/knowledge/show/lfunction.degree
[primitive]: https://beta.lmfdb.org/knowledge/show/lfunction.primitive
[conductor]: https://beta.lmfdb.org/knowledge/show/lfunction.conductor
[motivic weight]: https://beta.lmfdb.org/knowledge/show/lfunction.motivic_weight
[analytic conductor]: https://beta.lmfdb.org/knowledge/show/lfunction.analytic_conductor
[root analytic conductor]: https://beta.lmfdb.org/knowledge/show/lfunction.roo_analytic_conductor
[functional equation]: https://beta.lmfdb.org/knowledge/show/lfunction.functional_equation
[euler factors]: https://beta.lmfdb.org/knowledge/show/lfunction.euler_product
[arithmetic]: https://beta.lmfdb.org/knowledge/show/lfunction.arithmetic


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

