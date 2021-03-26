
Search table: lfunc_search
--------------------------

| Field | type | description |
|----------|    ------     | ----- |
origin                        |text                | **to be removed**
primitive                     |boolean             | true if L-func is [primitive], we use the second moment in many instances to decide this
conductor                     |numeric             | the [conductor] of the L-func
central_character             |text                | the conrey label of the primitive character that induces the [central character] of modulus equal the conductor
self_dual                     |boolean             | true if L-func is self-dual, i.e., all Dirichlet coefficients are real
motivic_weight                |smallint            | the [motivic weight] of the L-func
degree                        |smallint            | the [degree] of the L-func
order_of_vanishing            |smallint            | the [analytic rank], the order of vanishing at its central point
algebraic                     |boolean             | if the L-func is [arithmetic], i.e. normalized [Dirichlet coefficients] are algebraic numbers, conjecturally this is the same as being algebraic
z1                            |numeric             | the [lowest zero]
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
euler_factors                 |numeric[]           | **to be removed??**
rational                      |boolean             | if the [Dirichlet coefficients] in arithmetic normalisation are rational
euler{p}                      |numeric[]           | arrays of length degree + 1 representing the [euler factors] for p = 2, 3, 5,..., 97
root_analytic_conductor       |double precision    | the [root analytic conductor]
instance_types                |text[]              | representing the keys of the multimap url(type) -> url(instance)
instance_urls                 |text[]              | representing the values of the multimap url(type) -> url(instance)
spectral_label                |text                | the [spectral label], e.g. r0e3-p4.23p33.33m37.56
is_instance_{type}            |boolean             | type in {Artin,BMF,CMF,DIR,ECQ,ECQSymPower,ECNF,G2Q,HMF,MaassGL3,MaassGL4,MaassGSp4,NF}






Data table: lfunc_table
-----------------------
____________
| Field | type | description |
|----------|    ------     | ----- |
algebraic                     |boolean             | ^
analytic_conductor            |double precision    | ^
bad_lfactors                  |jsonb               | euler factors for the bad primes as an array of [p, [1, -ap, ...]]
bad_primes                    |bigint[]            | ^
central_character             |text                | ^
conductor                     |numeric             | ^
conductor_radical             |bigint              | ^
conjugate                     |text                | the [label] of the conjugate L-function, if self dual then None
degree                        |smallint            | ^
euler_factors                 |jsonb               | the first [euler factors] stored as array of lists, if the L-func is not rational, we represent each coefficient as pair of doubles corresponding to the real, imaginary pair. We need jsonb to be able to handle
euler_factors_factorization   |jsonb               | if the L-func is rational of degree larger than 4, we store the factorisation of the euler_factors
factors                       |text[]              | an array with the labels of the primitive factors (potentially repeated), where the last entry will be Null if we don't know the full factorization
index                         |smallint            | the last component of the label
instance_types                |text[]              | ^
instance_urls                 |text[]              | ^
label                         |text                | ^
leading_term_mid              |numeric             | the mid point of the leading term of the Taylor expansion of the L-function centered at t = 0 on the critical line
leaving_term_rad              |real                | the radious point of the leading term of the Taylor expansion of the L-function centered at t = 0 on the critical line
load_key                      |text                | a string marking the upload
motivic_weight                |smallint            | ^
mu_imag                       |numeric[]           | ^
mu_real                       |smallint[]          | ^
nu_imag                       |numeric[]           | ^
nu_real_doubled               |smallint[]          | ^
order_of_vanishing            |smallint            | ^
origin                        |text                | url for the object that was use to generated this data
plot_delta                    |float4              | the spacing of the plot_values
plot_values                   |float4[]            | the values of the Z function spaced by plot_delta, i.e., plot_values = [Z(k*plot_delta) for k in range(len(plot_delta))]
positive_zeros_mid            |numeric[]           | the midpoint of the first zeros, at most 10, represented as a ball
positive_zeros_rad            |double precision[]  | the radious of the first zeros, at most 10, represented as a ball
positive_zeros_extra          |double precision[]  | the remaining zeros via their correct double approximation
prelabel                      |text                | ^
primitive                     |boolean             | ^
rational                      |boolean             | ^
root_analytic_conductor       |double precision    | ^
root_angle                    |double precision    | ^
self_dual                     |boolean             | ^
special_values_at             |double precision[]  | an array of t's where L(t) has been computed as a ball, t is in analytic normalization
special_values_mid            |numeric[]           | the mid point of L(t)
special_values_rad            |real[]              | the radious of L(t)
spectral_label                |text                | ^
root_angle_mid                |numeric             | the mid point of the argument of [root number] normalized between -.5 to .5
root_angle_rad                |real                | the radious point of the argument of [root number] normalized between -.5 to .5
trace_hash                    |bigint              | ^
z1                            |numeric             | the first zero where all the last digit may have an error of +-1, e.g., we could represent pi as 3.1416
poles                         |double precision[]  | location of the poles in arithmetic normalization



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


TODO: add column is_instance_ECQSymPower
