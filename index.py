HEADER = "id|origin|primitive|conductor|central_character|self_dual|motivic_weight|Lhash|degree|order_of_vanishing|algebraic|z1|gamma_factors|trace_hash|root_angle|prelabel|analytic_conductor|mu_real|mu_imag|nu_real_doubled|nu_imag|bad_primes".split("|")
TYPES = "bigint|text|boolean|numeric|text|boolean|smallint|text|smallint|smallint|boolean|numeric|jsonb|bigint|double precision|text|double precision|smallint[]|numeric[]|smallint[]|numeric[]|bigint[]".split("|")

import re, sys
from six import string_types
from collections import Counter
from sage.all import psi, cached_function, ZZ, RR, GCD, ceil, RealField, ComplexField, CDF, walltime
try:
    from sage.rings.complex_mpfr import ComplexNumber
except ModuleNotFoundError:
    from sage.rings.complex_number import ComplexNumber
from sage.rings.real_mpfr import RealLiteral
#from lmfdb.backend.encoding import LmfdbRealLiteral
from dirichlet_conrey import DirichletGroup_conrey, DirichletCharacter_conrey

# This code comes from lmfdb.backend.encoding, but we don't want to import the lmfdb (in order to avoid db connection time)
class LmfdbRealLiteral(RealLiteral):
    """
    A real number that prints using the string used to construct it.
    """

    def __init__(self, parent, x=0, base=10):
        if not isinstance(x, string_types):
            x = str(x)
        RealLiteral.__init__(self, parent, x, base)

    def __repr__(self):
        return self.literal
def load(x, H, T):
    if x == r'\N':
        return None
    elif T == "text":
        return x
    elif T == "boolean":
        return True if x == "t" else False
    elif T in ["bigint", "smallint"] or H == "conductor":
        return ZZ(x)
    elif T in ["bigint[]", "smallint[]"]:
        return [ZZ(a) for a in x[1:-1].split(",") if a]
    elif T == "numeric[]":
        return [LmfdbRealLiteral(RR, a) for a in x[1:-1].split(",") if a]
    elif T == "double precision" or H == "z1":
        # Use LmfdbRealLiteral so that we can get the original string back
        return LmfdbRealLiteral(RR, x)
    elif H == "gamma_factors":
        # we don't care about gamma factors for the index
        return None
    else:
        raise RuntimeError((x, H, T))

def process_line(line, HEADER, TYPES):
    line = line.strip()
    L = {H: load(x, H, T) for (x, H, T) in zip(line.split("|"), HEADER, TYPES)}
    L['line'] = line
    return L


def find_prec(s):
    if isinstance(s, string_types):
        # strip negatives and exponent
        s = s.replace("-","")
        if "e" in s:
            s = s[:s.find("e")]
        return ceil(len(s) * 3.322)
    else:
        try:
            return s.parent().precision()
        except Exception:
            return 53


def hash_compare(x,y):
    if x['Lhash'] == y['Lhash']:
        return True
    if x['trace_hash'] is not None and y['trace_hash'] is not None:
        if x['trace_hash'] == y['trace_hash']:
            assert invariants_compare(x,y), "%s\n%s\n" % (x, y)
            return True
        else:
            return False
    return None

def invariants_compare(x, y):
    if x['primitive'] != y['primitive']:
        return False
    # check if z1s match
    z1s = [repr(elt['z1']) for elt in [x,y]]
    l = min([len(str(elt)) for elt in z1s])
    # we ignore the last 3 digits
    for k in range(l -3):
        if z1s[0][k] != z1s[1][k]:
            return False
    if x['order_of_vanishing'] != y['order_of_vanishing']:
        return False
    if (x['root_angle'] - y['root_angle']).abs() > 1e-7:
        return False
    # all the invariants agree
    return True

def compare(x, y):
    c = hash_compare(x, y)
    if c is None: #the Lhash differs, but might still be the same L-function
        if invariants_compare(x, y):
            return None
        else:
            return False
    else:
        return c


def compute_index(res):
    if len(res) == 1:
        res[0]['index'] = 0
    else:
        res.sort(key=lambda elt: elt['z1'])
        index = 0
        for i, x in enumerate(res):
            if 'index' in x:
                # already markes as unknown
                assert x['index'] == -1
                continue
            for j, y in enumerate(res):
                if i == j:
                    continue
                c = compare(x, y)
                if c is None:
                    x['index'] = -1
                    y['index'] = -1
                    break
                elif c:
                    if j > i:
                        # y doesn't have an index
                        continue
                    else: # J < i
                        assert invariants_compare(x, y)
                        print(x['origin'], y['origin'])
                        x['index'] = y['index']
                        break
            if 'index' not in x:
                x['index'] = index
                index += 1


def run(inputfilename, outputfilename):
    start = walltime()
    previouslabel = None
    res = []
    with open(inputfilename) as F:
        with open(outputfilename, 'w') as W:
            for i, line in enumerate(F.readlines(),1):
                if i == 0 and line.startswith('id|'):
                    HEADER = line.strip().split('|')
                    continue
                if i == 2 and line.startswith('bigint|'):
                    TYPES = line.strip().split('|')
                if i == 2 and not line:
                    continue

                L = process_line(line, HEADER, TYPES)
                L['line'] = line.strip()
                if L['prelabel'] == previouslabel:
                    res.append(L)
                else:
                    compute_index(res)
                    for r in res:
                        W.write('%s|%d\n' % (r['line'], r['index']))
                    res = [L]
                    previouslabel = L['prelabel']
                if i % 24200 == 0:
                    sys.stdout.write("%d (%.2f%%) lines done in %.2fs\tETA:%.2fs\n" % (i, i*100./24201375, walltime(start), (24201375 - i)*walltime(start)/i))
                    sys.stdout.flush()
            else:
                compute_index(res)
                for r in res:
                    W.write('%s|%d\n' % (r['line'], r['index']))


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])

