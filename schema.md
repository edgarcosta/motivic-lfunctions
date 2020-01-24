


Want to have a unified table, largely cnosisting exactly of what pops out of
Andy's code.




                            | Field      | type                                              | description |
                            | ---------- | ------                                            | -----       |
label                       | string     | ???
origin                      | string     | an LMFDB url from which we have derived this data
conductor                   | numeric
motivic_weight              | smallint
degree                      | smallint
mus                         | double[]   | the shifts of gamma_r
primitive                   | bool       |
root_number_arg             | double[]   | in some range
primitive_central_character | string     |
central_character           | string     | might be Null if not computed
self_dual                   | bool       |
conjugate_label             | string
load_key                    | string     | author/uploader ??
Lhash                       | string     | hash of lowest zero of L-function


ignore:
symmetry_type
