# Dedekin zeta function and artin Representations

We want the L-function of every Artin representation that may arise as the factor of the Dedekin zeta function of a field that has absolute discriminant less than < 10^6 or <10^10 if it is also 7-smooth, we expect about 4 to 5 million L-functions.

```
def count_artinrep_lfun(query):
    def chunkify(l, size=1000):
        res = [list(l)[elt*size:(elt + 1)*size] for elt in range(len(l)//size + 1)]
        return res
    res = {}
    missing = set({})
    for chunk in chunkify(list(db.nf_fields.search(query, 'coeffs'))):
        set_chunk = set(tuple(int(c) for c in elt) for elt in chunk)
        done = []
        for ar in db.artin_field_data.search({'$or': [{'Polynomial': elt} for elt in chunk]},
                                              ['ArtinReps', 'ConjClasses', 'Polynomial']):
            done.append(tuple(int(elt) for elt in ar['Polynomial']))
            for rep, data in zip(ar['ArtinReps'], ar['ConjClasses']):
                if int(data['Size']) != 1:
                    res[rep['Baselabel']] = int(data['Size'])
        missing = missing.union(set_chunk.difference(set(done)))
    return res, missing
res, miss = count_artinrep_lfun({'disc_abs' : {'$lte': 10^6}, 'degree': {'$gt': 2, '$lte': 20}})
print(sum(res.values()), len(miss))
3164727 223839
res, miss = count_artinrep_lfun({'disc_abs' : {'$lte': 10^6}, 'degree': {'$gt': 2, '$lte': 20}})
510156 19318
```



# Dirichlet Characters

We are aiming at
(conductor <= 1000 or (conductor <= 10000 and phi(order) <= 47)) ~ 4 million

Where the second clause guarantees that we get every Dirichlet character showing up as a factor in a desired Dedekin zeta function.

# Elliptic Curves

For every elliptic curve in the database we will have the L-function of Sym^d E if d <= 8 and Conductor(Sym^d E) <= 1e9.
This gives about ~ 4 million L-functions.

```
N <= 10^9 total: 4113861
     d    total
     1    2917287
     2    1073272
     3    40286
     4    64067
     5    655
     6    13389
     7    64
     8    4841
```


# Classical Modular forms

To fill in later.
