
Search table: lfunc_search
--------------------------

| Field | type | description |
|----------|    ------     | ----- |
primitive                     |boolean             | true if L-func is [primitive], we use the second moment in many instances to decide this
conductor                     |numeric             | the [conductor] of the L-func
central_character_prim        |text                | the conrey label of the primitive character that induces the [central character] of modulus equal the conductor
self_dual                     |boolean             | true if L-func is self-dual, i.e., all Dirichlet coefficients are real
motivic_weight                |smallint            | the [motivic weight] of the L-func
degree                        |smallint            | the [degree] of the L-func
order_of_vanishing            |smallint            | the [analytic rank], the order of vanishing at its central point
algebraic                     |boolean             | if the L-func is [arithmetic], i.e. normalized [Dirichlet coefficients] are algebraic numbers, conjecturally this is the same as being algebraic
z1                            |numeric             | the [lowest zero] (the last digit may have an error of +-1, e.g., we could represent pi as 3.1416)
trace_hash                    |bigint              | linear combination of the a_p between 2^12 and 2^13 reduced mod 2^61-1 as defined in Section 4.3 of [BSSVY](https://arxiv.org/abs/1602.03715), only for rational L-functions
root_angle                    |double precision    | the argument of [root number] normalized between -.5 to .5
prelabel                      |text                | the label without the index
analytic_conductor            |double precision    | the [analytic conductor]
mu_real                       |smallint[]          | the real part (in [0, 1]) of mus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
mu_imag                       |numeric[]           | the imaginary part of mus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
nu_real_doubled               |smallint[]          | the real part of nus doubled, so they are integers, in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
nu_imag                       |numeric[]           | the imaginary part of nus in the analytic normalization the [functional equation], where if possible Gamma_R factors have been converted to Gamma_C
bad_primes                    |bigint[]            | primes dividing the [conductor]
label                         |text                | the [label] of the lfunc, e.g.: 3-1-1.1-r0e3-p4.23p33.33m37.56-0
index                         |smallint            | the last component of the [label]
conductor_radical             |bigint              | product of the bad_primes, i.e., the primes dividing the [conductor]
dirichlet_coefficients        |numeric[]           | the [Dirichlet coefficients] for rational L-functions in arithmetic normalisation starting with a_1
rational                      |boolean             | if the [Dirichlet coefficients] in arithmetic normalisation are rational
euler{p}                      |numeric[]           | arrays of length degree + 1 representing the [euler factors] for p = 2, 3, 5,..., 97 (only for rational L-functions)
root_analytic_conductor       |double precision    | the [root analytic conductor]
instance_types                |text[]              | representing the keys of the multimap url(type) -> url(instance)
instance_urls                 |text[]              | representing the values of the multimap url(type) -> url(instance)
spectral_label                |text                | the [spectral label], e.g. r0e3-p4.23p33.33m37.56
is_instance_{type}            |boolean             | type in {Artin,BMF,CMF,DIR,ECQ,ECQSymPower,ECNF,G2Q,HMF,MaassGL3,MaassGL4,MaassGSp4,NF,HGM}






Data table: lfunc_table
-----------------------
____________
| Field | type | description |
|----------|    ------     | ----- |
bad_lfactors                  |jsonb               | Euler factors for the bad primes as an array of [p, [1, -ap, ...]], if p is larger than 100
conjugate                     |text                | the [label] of the conjugate L-function, if self-dual, then None
euler_factors                 |jsonb               | the first [euler factors] stored as array of lists, if the L-func is not rational, we represent each coefficient as pair of doubles corresponding to the real, imaginary pair. We need jsonb to support different formats in the same column.  Store for p < 100.
euler_factors_factorization |jsonb | If the L-func is rational of a degree larger than 4, we store the factorization of the euler_factors.  Need jsonb since factors can have different lengths
factors                       |text[]              | an array with the labels of the primitive factors (potentially repeated), where the last entry will be Null if we don't know the full factorization
factors_shift                 |numeric[]           | store the array s where L(s) = L1(s+w1) L2(s + w2) ... Ln(s + wn), NULL if all zeros
index                         |smallint            | the last component of the label
label                         |text                | the [label] of the lfunc, e.g.: 3-1-1.1-r0e3-p4.23p33.33m37.56-0
leading_term_mid              |numeric             | the midpoint of the leading term of the Taylor expansion of the L-function centered at t = 0 on the critical line
leading_term_rad              |float8              | the radius point of the leading term of the Taylor expansion of the L-function centered at t = 0 on the critical line
origin                        |text                | URL for the object that was used to generate this data
plot_x                        |float4[]            | x-coordinates for local minima and local maxima of the [Z-function] (only for non-negative x)
plot_y                        |float4[]            | y-coordinates for the local minima and local maxima of the [Z-function] at the values of plot_x
plot_deriv                    |float4[]            | a list of length equal to the list of zeros given in positive_zeros_mid and positive_zeros_extra, giving the first derivative at each 0
plot_extra                    |float4[]            | a list of triples [x,y,z], constraining y=Z(x) and z=Z^(1)(x)
positive_zeros_mid            |numeric[]           | the midpoint of the first zeros, at most 10, represented as a ball
positive_zeros_rad            |float8[]            | the radius of the first zeros, at most 10, represented as a ball
positive_zeros_extra          |float8[]            | the remaining zeros via their correct double approximation.  These should be stored as far as we want the plot to extend.
root_angle_mid                |numeric             | the midpoint of the argument of [root number] normalized between -.5 to .5
root_angle_rad                |float8              | the radius point of the argument of [root number] normalized between -.5 to .5
special_values_at             |float8[]            | an array of {ti} where L^(ki)(ti) has been computed as a ball, ti is in arithmetic normalization
special_values_der_order      |smallint[]          | the ki on the line above
special_values_mid            |numeric[]           | the midpoint of L^(ki)(ti), stored either as a pair (giving real and imaginary part) or a single float (if all special values are real)
special_values_rad            |float8[]            | the radius of L^(ki)(ti)


[Dirichlet coefficient]: https://beta.lmfdb.org/knowledge/show/lfunction.dirichlet_series
[analytic conductor]: https://beta.lmfdb.org/knowledge/show/lfunction.analytic_conductor
[analytic rank]: https://beta.lmfdb.org/knowledge/show/lfunction.analytic_rank
[arithmetic]: https://beta.lmfdb.org/knowledge/show/lfunction.arithmetic
[central character]: https://beta.lmfdb.org/knowledge/show/lfunction.central_character
[conductor]: https://beta.lmfdb.org/knowledge/show/lfunction.conductor
[degree]: https://beta.lmfdb.org/knowledge/show/lfunction.degree
[euler factors]: https://beta.lmfdb.org/knowledge/show/lfunction.euler_product
[functional equation]: https://beta.lmfdb.org/knowledge/show/lfunction.functional_equation
[label]: https://beta.lmfdb.org/knowledge/show/lfunction.label
[lowest zero]: https://beta.lmfdb.org/knowledge/show/lfunction.zeros
[motivic weight]: https://beta.lmfdb.org/knowledge/show/lfunction.motivic_weight
[primitive]: https://beta.lmfdb.org/knowledge/show/lfunction.primitive
[root analytic conductor]: https://beta.lmfdb.org/knowledge/show/lfunction.root_analytic_conductor
[root number]: https://beta.lmfdb.org/knowledge/show/lfunction.sign
[spectral label]: https://beta.lmfdb.org/knowledge/show/lfunction.spectral_label
[Z-function]: https://beta.lmfdb.org/knowledge/show/lfunction.zfunction

Old plot info
-------------

* plot_delta (float4): the spacing of the plot_values
* plot_values (float4[]): the values of the Z function spaced by plot_delta, i.e., plot_values = [Z(k*plot_delta) for k in range(len(plot_delta))]

