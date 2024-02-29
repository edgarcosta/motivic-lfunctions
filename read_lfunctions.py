import ast
import fcntl
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from argparse import ArgumentParser
from collections import Counter, defaultdict, OrderedDict
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from dirichlet_conrey import DirichletGroup_conrey, DirichletCharacter_conrey
from halo import Halo
from itertools import islice
from sage.all import (
    CDF,
    ComplexBallField,
    EllipticCurve,
    GCD,
    GF,
    HyperellipticCurve,
    Integer,
    Integers,
    NumberField,
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    RBF,
    RDF,
    RR,
    RealBallField,
    RealIntervalField,
    ZZ,
    cached_function,
    ceil,
    gcd,
    gp,
    lazy_attribute,
    lazy_class_attribute,
    log,
    magma,
    next_prime,
    prime_powers,
    prime_range,
    primes_first_n,
    prod,
    psi,
    spline,
    srange,
    vector,
)
from sage.rings.real_arb import RealBall
from sage.rings.real_double import RealDoubleElement
from sage.rings.real_mpfr import RealLiteral, RealField, RealNumber
from sage.structure.unique_representation import CachedRepresentation


def mod1(elt):
    return RDF(elt) - RDF(elt).floor()


lhash_regex = re.compile(r'^\d+([,]\d+)*$')
ec_lmfdb_label_regex = re.compile(r'(\d+)\.([a-z]+)(\d*)')


def progress_bar(current, total, time, barLength=10):
    percent = float(current) * 100 / total
    if current > 0:
        eta = '%.2f' % (time * (total - current) / float(current),)
    else:
        eta = '+oo'
    done = '▰' * int(percent / 100 * barLength)
    notdone = '▱' * (barLength - len(done))

    return ' Progress: [%s%s] %.2f%% ETA: %s seconds' % (done, notdone, percent, eta)


def chunkify(l, size=100):
    return [islice(l, elt * size, (elt + 1) * size)
            for elt in range(len(l) // size + 1)]


def coeff_to_poly(c, var=None):
    # from lmfdb/utils/utilities.py
    """
    Convert a list or string representation of a polynomial to a sage polynomial.

    Examples:
    >>> coeff_to_poly("1 - 3x + x^2")
    x**2 - 3*x + 1
    >>> coeff_to_poly("1 - 3*x + x**2")
    x**2 - 3*x + 1
    """
    if isinstance(c, str):
        # accept latex
        c = c.replace("{", "").replace("}", "")
        # autodetect variable name
        if var is None:
            varposs = set(re.findall(r"[A-Za-z_]+", c))
            if len(varposs) == 1:
                var = varposs.pop()
            elif not(varposs):
                var = 'x'
            else:
                raise ValueError("Polynomial must be univariate")
    if var is None:
        var = 'x'
    return PolynomialRing(QQ, var)(c)


def force_list_int(elt):
    if isinstance(elt, Iterable):
        return [force_list_int(x) for x in elt]
    assert isinstance(elt, (int, Integer))
    return int(elt)


# READ utils
def strip(elt):
    return elt.replace(' ', '').strip()


def atoii(elt, level=1, stripped=True):
    s = strip(elt) if not stripped else elt
    assert level >= 1
    if set('[]') == set(s):
        return ast.literal_eval(s)
    assert s[:level] == '[' * level and s[-level:] == ']' * level, s
    if level == 1:
        return list(map(int, s[1:-1].split(',')))
    else:
        s = s[level:-1].split('[' * (level - 1))
        return [atoii('[' * (level - 1) + c.rstrip(','), level=level - 1) for c in s]


# WRITE utils
def json_dumps(inp, typ, recursing=False):
    if inp is None:
        return u"\\N"
    elif typ in ("text", "char", "varchar"):
        if not isinstance(inp, str):
            inp = str(inp)
        inp = (
            inp.replace("\\", "\\\\")
            .replace("\r", r"\r")
            .replace("\n", r"\n")
            .replace("\t", r"\t")
            .replace('"', r"\"")
        )
        if recursing and ("{" in inp or "}" in inp):
            inp = '"' + inp + '"'
        return inp
    elif typ in ("json", "jsonb"):
        # we are only using this for lists of lists of integers and None
        return json.dumps(inp)
    elif isinstance(inp, str):
        # assuming that we already took care of the right conversion
        # as typ != text, char, varchar
        return inp
    elif typ[-2:] == "[]":
        if not isinstance(inp, (list, tuple)):
            raise TypeError("You must use list or tuple for array columns")
        if not inp:
            return "{}"
        subtyp = None
        sublen = None
        for x in inp:
            if isinstance(x, (list, tuple)):
                if subtyp is None:
                    subtyp = typ
                elif subtyp != typ:
                    raise ValueError("Array dimensions must be uniform")
                if sublen is None:
                    sublen = len(x)
                elif sublen != len(x):
                    raise ValueError("Array dimensions must be uniform")
            elif subtyp is None:
                subtyp = typ[:-2]
            elif subtyp != typ[:-2]:
                raise ValueError("Array dimensions must be uniform")
        return "{" + ",".join(json_dumps(x, subtyp, recursing=True) for x in inp) + "}"
    elif isinstance(inp, (float, int, Integer, RealNumber, RealLiteral, RealDoubleElement)):
        return str(inp).replace("L", "")
    elif typ == "boolean":
        return "t" if inp else "f"
    else:
        raise TypeError("Invalid input %s (%s) for postgres type %s" % (inp, type(inp), typ))


float_re = r'[-+]?(?:[0-9]*[.]?[0-9]+(?:[ed][-+]?[0-9]+)?|inf|nan)'
FLOAT_RE = re.compile(float_re)
ARB_RE = re.compile(r'\[({float_re}) \+/- ({float_re})\]'.format(
    float_re=float_re))


def realball_to_mid_rad_str(elt):
    # we really on arb_printn
    # https://arblib.org/arb.html#c.arb_printn
    s = str(elt)
    if elt.rad() == 0:
        mid, rad = s, '0'
    else:
        match = ARB_RE.match(s)
        mid, rad = match[1], match[2]
    assert elt in elt.parent()(r'[{mid} +/- {rad}]'.format(mid=mid, rad=rad))
    return mid, rad


def complexball_to_mid_rad_str(elt, extra_digits=9):
    m1, r1 = realball_to_mid_rad_str(elt.real())
    m2, r2 = realball_to_mid_rad_str(elt.imag())
    return '{%s,%s}' % (m1, m2), '{%s,%s}' % (r1, r2)


class RealBallLiteral(RealBall):
    def __init__(self, parent, mid, rad):
        self.literal = r'[{mid} +/- {rad}]'.format(mid=mid, rad=rad)
        return parent(self.literal)

    def __repr__(self):
        return self.literal


def enough_digits(elt, digits=6):
    return elt.abs().rad_as_ball() * 10**digits - 1 < 0


def mid_point_100bits(elt):
    assert elt.rad().log(2) <= -100
    elt100 = elt * 2**100
    eltint = elt100.mid().round()
    ball = elt.parent()([eltint - 1, eltint + 1]) << -100
    assert elt in ball
    return eltint, ball


def ball_from_midpoint(elt):
    assert isinstance(elt, str)
    rn = RealField(200)(elt)
    eltint = (rn << 100).round()
    ball = RealBallField(200)([eltint - 1, eltint + 1]) << -100
    return ball


LOG_TEN_TWO_PLUS_EPSILON = float.fromhex('0x3.5269e12f346e4p+0')


def numeric_to_ball(elt):
    """
    converts a string representing a real number into a ball
    by setting the radius such that the last two digits are unknown
    as one usually doesn't truncate when printing floats to be able to recover all the binary digits

    sage: RIF(numeric_to_ball('1.2500'))
    1.25?
    sage: numeric_to_ball('1.2500').endpoints()
    (1.24199999998381, 1.25800000001619)
    """
    if isinstance(elt, float) or elt in RDF:
        elt = RR(elt)
        return RBF(elt, elt.ulp() * 0.5)
    assert isinstance(elt, str)
    sigfig_mantissa = elt.lstrip('-0.')
    sigfigs = len(sigfig_mantissa) - ('.' in sigfig_mantissa)
    bits = int(LOG_TEN_TWO_PLUS_EPSILON * sigfigs) + 1
    # note that the precision is already higher than what elt represents
    assert '.' in elt
    # the last 3 digits might be off
    rad = 10**(-(len(elt.split('.')[-1]) - 3))
    return RealBallField(bits)(elt, rad)


def realnumber_to_ball(elt, R):
    return R(elt, float(elt.ulp()))


def approx_ball(elt, prec=53):
    """
    if we can approximate the ball, returns such approximation
    """
    # this is what we would get from such approximation
    approx_ball = realnumber_to_ball(elt.numerical_approx(prec=prec), RealBallField(prec))
    if elt in approx_ball:
        return approx_ball.mid()
    else:
        return None


# to avoid the discontinuity at (-inf, 0], which will result in
# [+/- 3.15]
def arg_hack(foo):
    if not foo.real().contains_zero() and foo.real().mid() < 0:
        arg = (-foo).arg()
        if arg > 0:
            arg -= foo.parent().pi().real()
        else:
            arg += foo.parent().pi().real()
        return arg
    else:
        return foo.arg()


def normalized_arg(foo):
    arg = arg_hack(foo) / (2 * foo.parent().pi())
    while arg > 0.5:
        arg -= 1
    while arg <= -0.5:
        arg += 1
    return arg


def extend_multiplicatively(Z):
    for pp in prime_powers(len(Z) - 1):
        for k in range(1, (len(Z) - 1) // pp + 1):
            if gcd(k, pp) == 1:
                Z[pp * k] = Z[pp] * Z[k]


def dirichlet_coefficients(euler_factors):
    R = vector(sum(euler_factors), []).base_ring()
    PS = PowerSeriesRing(R)
    pef = list(zip(primes_first_n(len(euler_factors)), euler_factors))
    an_list_bound = next_prime(pef[-1][0])
    res = [1] * an_list_bound
    for p, ef in pef:
        k = RR(an_list_bound).log(p).floor() + 1
        foo = (1 / PS(ef)).padded_list(k)
        for i in range(1, k):
            res[p**i] = foo[i]
    extend_multiplicatively(res)
    return res


@cached_function
def DirGroup(m):
    return DirichletGroup_conrey(m)


@cached_function
def primitivize(label):
    m, n = [ZZ(a) for a in label.split(".")]
    char = DirichletCharacter_conrey(DirGroup(m), n).primitive_character()
    return "%d.%d" % (char.modulus(), char.number())


def prod_central_character(labels):
    char = prod([
        DirichletCharacter_conrey(DirGroup(m), n).primitive_character()
        for m, n in [[ZZ(a) for a in label.split(".")] for label in labels]
    ]).primitive_character()
    return "%d.%d" % (char.modulus(), char.number())


logpi = RR.pi().log()
log2pi = (2 * RR.pi()).log()


def log_L_inf(s, mu, nu):
    return ((sum([psi((s + elt) / 2) for elt in mu]) -
             len(mu) * logpi) / 2 +
            sum([psi(s + elt) for elt in nu]) -
            len(nu) * log2pi)


@cached_function
def conductor_an(GR, GC):
    return (2 * log_L_inf(1 / 2, GR, GC).real()).exp()


class LfunctionsParser(object):
    default_prec = 300
    sep = ':'

    @lazy_class_attribute
    def CBF(cls):
        return ComplexBallField(cls.default_prec)

    @lazy_class_attribute
    def RBF(cls):
        return RealBallField(cls.default_prec)

    def __init__(self, in_headers, out_headers):
        self.out_headers = out_headers
        self.in_headers = in_headers

    def from_arb2(self, lower, upper, exp):
        return self.RBF([lower, upper]) << exp

    def from_acb2(self, lower_real, upper_real, exp_real, lower_imag, upper_imag, exp_imag):
        return self.CBF(self.RBF([lower_real, upper_real]) << exp_real, self.RBF([lower_imag, upper_imag]) << exp_imag)

    def read_vector_float(self, elt):
        if elt == '[]':
            return []
        return [RDF(x) for x in elt[1:-1].split(',')]

    def read_origin_label(self, elt):
        return elt

    def read_root_number_acb(self, elt):
        return self.from_acb2(*sum(atoii(elt, level=2), []))

    def read_order_of_vanishing(self, elt):
        return int(elt)

    def read_leading_term_arb(self, elt):
        return self.from_arb2(*atoii(elt))

    def read_special_values_acb(self, elt):
        values = [self.from_acb2(*sum(x, [])) for x in atoii(elt, level=3)]
        return dict(enumerate(values, 1))

    def read_positive_zeros_arb(self, elt):
        return [self.from_arb2(*x) for x in atoii(elt, level=2)]

    def read_positive_zeros_extra(self, elt):
        return self.read_vector_float(elt)

    def read_plot_delta(self, elt):
        return RDF(elt)

    def read_plot_values(self, elt):
        return self.read_vector_float(elt)

    def read_mus(self, elt):
        return self.read_vector_float(elt)

    def read_conductor(self, elt):
        return int(elt)

    def read_trace_hash(self, elt):
        return int(elt)

    def read_second_moment(self, elt):
        return RDF(elt).round()

    def read_line(self, elt, headers):
        return {k: getattr(self, 'read_' + k)(v) for k, v in zip(headers, strip(elt).split(self.sep))}

    def read_out_line(self, line):
        return self.read_line(line, self.out_headers)

    def read_in_line(self, line):
        return self.read_line(line, self.in_headers)

    def read_record(self, in_line, out_line):
        raise NotImplementedError("you should override this ;)")

    def read_files(self, in_file, out_file, begin=0, end=-1):
        print(in_file, out_file, begin, end)
        with open(in_file) as iF:
            with open(out_file) as oF:
                for i, (in_line, out_line) in enumerate(zip(iF, oF)):
                    if i < begin:
                        continue
                    if i == end:
                        break
                    for elt in self.process_record(self.read_record(in_line, out_line)):
                        yield elt

    def process_record(self, rec):
        yield rec


class lfunction_element(object):
    # desired columns when creating object from db
    projection = [
        'Lhash',
        'accuracy',
        'central_character',
        'conductor',
        'degree',
        'id',
        'index',
        'label',
        'motivic_weight',
        'order_of_vanishing',
        'positive_zeros',
        'positive_zeros_mid',
        'positive_zeros_rad',
        'prelabel',
        'primitive',
        'root_analytic_conductor',
        'root_angle',
        'trace_hash',
    ]
    euler_factors_bound = 100
    RBF = LfunctionsParser.RBF
    CBF = LfunctionsParser.CBF
    ZZT = PolynomialRing(ZZ, "T")

    def __init__(self, data, from_db=False):
        self.algebraic = True
        self.coeff_info = None
        self.dirichlet_coefficients = None  # we prefer ai and euler_factors
        self.credit = None
        self.group = None
        self.st_group = None  # we don't have the infrastructure to determine this on the fly
        self.symmetry_type = None
        self.sign_arg = None
        if from_db:
            # convert from literal to integer
            data['conductor'] = int(data['conductor'])
            if data.get('instance_types') is None:
                data['instance_types'] = []
            if data.get('instance_urls') is None:
                data['instance_urls'] = []
        self.from_db = from_db

        self.__dict__.update(data)

        if 'leading_term_arb' in data:
            self.leading_term_mid, self.leading_term_rad = realball_to_mid_rad_str(self.leading_term_arb)
            self.leading_term = None

        if 'root_number_acb' in data:
            # sets
            # - root_angle: float
            # - sign_arg_mid: numeric
            # - sign_arg_rad: real
            # - root_number_mid: numeric[]
            # - root_number_rad: numeric[]
            self.root_number = None
            self.sign_arg = None
            if self.self_dual:
                # pindown root_number and its argument
                assert self.root_number_acb.contains_integer()
                root_number = Integer(self.root_number_acb)
                # pin it down
                self.root_number_acb = self.root_number_acb.parent()(root_number)
                arg = self.root_number_acb.real().parent()(0 if root_number == 1 else 0.5)
            else:
                arg = normalized_arg(self.root_number_acb)

            # self.root_number_mid,  self.root_number_rad = complexball_to_mid_rad_str(self.root_number_acb)
            self.root_angle_mid, self.root_angle_rad = realball_to_mid_rad_str(arg)

        if 'special_values_acb' in data:
            if self.self_dual:
                # pin down the values to the real axis
                for t, elt in self.special_values_acb.items():
                    assert elt.imag().contains_zero()
                    self.special_values_acb[t] = elt.parent()(elt.real())
            sv = [(t, complexball_to_mid_rad_str(elt))
                  for (t, elt) in sorted(self.special_values_acb.items()) if enough_digits(elt)]
            # at the moment artin and smalljac compute at the same points
            self.special_values_at = '{%s}' % ','.join(str(elt[0])for elt in sv)
            self.special_values_mid = '{%s}' % ','.join(elt[1][0] for elt in sv)
            self.special_values_rad = '{%s}' % ','.join(elt[1][1] for elt in sv)
            self.values = None

        if 'positive_zeros_arb' in data:
            self.accuracy = self.precision = None  # we are using balls
            # we increase the radius to 2**-103 => 31 digits
            assert(len(self.positive_zeros_arb) <= 10)
            max_rad = 2**-103
            # doesn't change the radius if we pass a negative value
            z = [realball_to_mid_rad_str(elt.add_error(max_rad - elt.rad()))
                 for elt in self.positive_zeros_arb]
            self.positive_zeros_mid = '{%s}' % ','.join([elt[0] for elt in z])
            self.positive_zeros_rad = '{%s}' % ','.join([elt[1] for elt in z])
            R = RealIntervalField(self.positive_zeros_arb[0].mid().prec())
            # remove the question mark and the unknown digit
            self.z1 = str(R(self.positive_zeros_arb[0]))[:-2]
            self.z2 = str(R(self.positive_zeros_arb[1]))[:-2]
            self.z3 = str(R(self.positive_zeros_arb[2]))[:-2]

    @lazy_attribute
    def positive_zeros_arb(self):
        if self.from_db:
            # we will use zero balls for comparisons to figure out label and/or factors
            if getattr(self, 'positive_zeros_rad', None) and getattr(self, 'positive_zeros_mid', None):
                R = RealBallField(max(53, self.positive_zeros_mid[0].prec() + 10))
                return [RealBallLiteral(R, m, r) for m, r in zip(self.positive_zeros_mid, self.positive_zeros_rad)]
            elif getattr(self, 'positive_zeros', None):
                if getattr(self, 'accuracy', None) == 100:
                    positive_zeros_arb = [ball_from_midpoint(elt) for elt in self.positive_zeros]
                    assert self.Lhash.split(',')[0] == str(mid_point_100bits(positive_zeros_arb[0])[0])
                else:
                    positive_zeros_arb = [numeric_to_ball(elt) for elt in self.positive_zeros]
                self.positive_zeros_extra = []
                return positive_zeros_arb
            else:
                assert False

    @lazy_attribute
    def positive_zeros_extra(self):
        return []

    @lazy_attribute
    def positive_zeros_extra_arb(self):
        return [numeric_to_ball(elt) for elt in self.positive_zeros_extra]

    def __repr__(self):
        if getattr(self, 'label', None):
            return "L-function %s" % self.label
        if hasattr(self, 'origin_label'):
            return "L-function obtained from %s" % self.origin_label
        return "L-function with degree %d and conductor %d" % (self.degree, self.conductor)

    def __eq__(self, other, pedantic=True):
        # first compare exact invariants
        for attr in ['prelabel', 'primitive', 'order_of_vanishing', 'motivic_weight']:
            if getattr(self, attr) != getattr(other, attr):
                return False

        # before comparing zeros check if type and origins match
        if self.origin == other.origin and self.origin:
            # the easy invariants check passed so if the origins are equal
            # the trace hash should be sufficient as comparison
            if self.trace_hash and other.trace_hash:
                return self.trace_hash == other.trace_hash

        # non exact invariants
        if mod1(self.root_angle - other.root_angle) > 1e-7:
            return False

        # check if the first zeros overlap
        # we could just check the first couple, but the dominant cost is the instantiation of the balls
        for i, (z, w) in enumerate(zip(self.positive_zeros_arb, other.positive_zeros_arb)):
            if not z.overlaps(w):
                if self.primitive and i > 0:
                    assert False, "we matched the first zero but not the %d-th" % (i + 1)
                return False

        # compare hashes

        # no every L-function uses the Lhash based on z1 of its factors
        if self.Lhash is not None and other.Lhash is not None and lhash_regex.match(self.Lhash) and lhash_regex.match(other.Lhash):
            # the lhash_regex also catches Lhash based on trace_hash
            # some we might be comparing trace_hash vs hash based on z1
            if self.Lhash == other.Lhash:
                return True

        if self.trace_hash is not None and other.trace_hash is not None:
            if self.trace_hash != other.trace_hash:
                return False
            # the trace_hashes match, but we might have a collision due to high degree
            if self.degree <= 20:
                # the lowest degree collision we are aware is at degree 32
                # https://beta.lmfdb.org/ModularForm/GL2/Q/holomorphic/3381/1/o/c/
                # https://beta.lmfdb.org/ModularForm/GL2/Q/holomorphic/3381/1/c/b/
                return True

        if pedantic:
            raise ValueError("Couldn't decide if the L-functions are equal %s, %s" % (self.__dict__, other.__dict__))
        else:
            return True

    def compute_factorization(self, possible_factors):
        # we are assuming that it factors with no repeated factors
        if not possible_factors:
            raise ValueError("possible_factors cannot be empty")

        if not isinstance(possible_factors, list):
            raise ValueError("possible_factors must be a list")

        # quotient properties

        conductor = self.conductor
        degree = self.degree
        order_of_vanishing = self.order_of_vanishing
        zeros = self.positive_zeros_arb + self.positive_zeros_extra_arb

        def divides(elt):
            if not (conductor % elt.conductor == 0 and
                    elt.motivic_weight <= self.motivic_weight and
                    elt.degree <= degree and
                    elt.order_of_vanishing <= order_of_vanishing):
                return False

            # check if the first 3 zeros show up in other
            k = 0
            for j, z in enumerate(elt.positive_zeros_arb[:3]):
                for i, w in enumerate(zeros[k:], k):
                    # all elements of w are less than all the elements of z
                    if w < z:
                        continue
                    # z overlaps with w, we can move to the next zero
                    if w.overlaps(z):
                        k = i + 1
                        break
                    # all elements of w are greater than all the elements of z
                    # and thus it will also be true for every other future w
                    if w > z:
                        if j > 0:
                            assert False, "we matched the first zero but not the %d-th" % (j + 1)
                        return False
                else:
                    # we got to the end of zeros
                    # and they were all smaller than z
                    assert j > 0  # assuring at least a zero matched
                    break
            return True

        def divide(elt):
            new_conductor = conductor // elt.conductor
            new_degree = degree - elt.degree
            new_order_of_vanishing = order_of_vanishing - elt.order_of_vanishing
            new_zeros = []
            elt_zeros = elt.positive_zeros_arb + elt.positive_zeros_extra_arb
            elt_zeros.reverse()
            k = 0
            while elt_zeros:
                w = elt_zeros.pop()
                for i, z in enumerate(zeros[k:], k):
                    # all elements of w are greater than all the elements of z
                    if z < w:
                        new_zeros.append(z)

                    # z overlaps with w, we can move to the next zero
                    if z.overlaps(w):
                        k = i + 1
                        break
                    if z > w:
                        assert False, 'we could not find the zero %s in %s' % (w, zeros[k:])
                else:
                    # got to the end and every element is larger than the current w
                    # and therefore we can't remove any more zeros
                    elt_zeros = []
            assert len(new_zeros) < len(zeros)
            return new_degree, new_conductor, new_order_of_vanishing, new_zeros

        res = []
        if isinstance(possible_factors[0], list):
            # we are given a possible list of factors for each factor
            for factors in possible_factors:
                for elt in factors:
                    elt = lfunction_element(elt, from_db=True)
                    if divides(elt):
                        degree, conductor, order_of_vanishing, zeros = divide(elt)
                        res.append(elt)
                        break
        else:
            for elt in possible_factors:
                elt = lfunction_element(elt, from_db=True)
                if divides(elt):
                    degree, conductor, order_of_vanishing, zeros = divide(elt)
                    res.append(elt)
                    if degree == 0:
                        break

        res = sorted(res, key=lambda x: (x.degree, x.conductor, x.positive_zeros_arb[0]))
        if degree == 0:
            self.factors_obj = res
            for elt in ['factors', 'Lhash', 'LhashArray']:
                if hasattr(self, elt):
                    self.__dict__.pop(elt)
        else:
            self.Lhash = "_" + str(mid_point_100bits(self.positive_zeros_arb[0])[0])
            self.LhashArray = [self.Lhash]
            self.factors = [elt.label for elt in res] + [None]

    def merge_instances(self, other):
        # merge instances
        assert len(other.instance_types) == len(other.instance_urls)
        m = set(zip(self.instance_types + other.instance_types,
                    self.instance_urls + other.instance_urls))
        self.instance_types, self.instance_urls = zip(*m)
        other.instance_types = self.instance_types
        other.instance_urls = self.instance_urls

    @classmethod
    def from_factors(cls, factors, data={}, **kwargs):
        """
        one may pass precomputed data via optional data parameter, e.g., euler factors
        """
        assert len(factors) >= 1
        if len(factors) == 1:
            return factors[0]
        factors.sort(key=lambda elt: (elt.degree, elt.conductor, elt.positive_zeros_arb[0]))
        data['degree'] = sum(elt.degree for elt in factors)
        data['conductor'] = prod(elt.conductor for elt in factors)
        data['factors_obj'] = factors
        data['motivic_weight'] = factors[0].motivic_weight
        assert all(data['motivic_weight'] == elt.motivic_weight for elt in factors)
        data['order_of_vanishing'] = sum(elt.order_of_vanishing for elt in factors)
        data['poles'] = sorted(set(sum([elt.poles for elt in factors], [])))
        data['primitive'] = False

        # Handle the zeros
        positive_zeros_arb = sum([elt.positive_zeros_arb for elt in factors], [])
        R = positive_zeros_arb[0].parent()
        positive_zeros_arb += sum([[realnumber_to_ball(z, R) for z in elt.positive_zeros_extra] for elt in factors], [])
        positive_zeros_arb.sort()
        data['positive_zeros_arb'] = positive_zeros_arb[:10]
        data['positive_zeros_extra'] = []
        rh_limit = 64 / data['degree']
        for elt in positive_zeros_arb[10:]:
            approx = approx_ball(elt)
            if elt is None or elt.mid() > rh_limit:
                break
            else:
                data['positive_zeros_extra'].append(approx)

        if all(getattr(elt, 'self_dual', False) for elt in factors):
            data['self_dual'] = True

        if all(getattr(elt, 'central_character', None) for elt in factors):
            data['central_character'] = prod_central_character([elt.central_character for elt in factors])

        if all(hasattr(elt, 'leading_term_arb') for elt in factors):
            data['leading_term_arb'] = prod(elt.leading_term_arb for elt in factors)

        if all(hasattr(elt, 'root_number_acb') for elt in factors):
            data['root_number_acb'] = prod(elt.root_number_acb for elt in factors)
        else:
            # we will use this for comparisons
            arg = mod1(sum(elt.root_angle for elt in factors))
            while arg > 0.5:
                arg -= 1
            while arg <= -0.5:
                arg += 1
            data['root_angle'] = arg

        if all(hasattr(elt, 'special_values_acb') for elt in factors):
            at = set(factors[0].special_values_acb)
            for elt in factors[1:]:
                at.intersection_update(elt.special_values_acb)
            data['special_values_acb'] = {t: prod(elt.special_values_acb[t] for elt in factors) for t in at}

        if all(hasattr(elt, 'plot_delta') and hasattr(elt, 'plot_values') for elt in factors):
            factor_plot_values = [[(j * elt.plot_delta, z) for j, z in enumerate(elt.plot_values)] for elt in factors]
            interpolations = [spline(elt) for elt in factor_plot_values]
            plot_range = 64.0 / data['degree']
            # we cannot hope to get a finer resolution
            data['plot_delta'] = plot_delta = max(plot_range / 256, max(elt.plot_delta for elt in factors))
            # we also don't want to extrapolate data
            assert all(plot_range <= elt[-1][0] for elt in factor_plot_values), '%s %s' % (plot_range, [elt[-1][0] for elt in factor_plot_values])
            data['plot_values'] = [prod([elt(i) for elt in interpolations]) for i in srange(0, plot_range + plot_delta, plot_delta)]
            assert len(data['plot_values']) <= 257

        if all(hasattr(elt, 'gamma_factors') for elt in factors):
            data['gamma_factors'] = [sorted(sum([elt.gamma_factors[i] for elt in factors], []))
                                     for i in range(2)]
        if all(hasattr(elt, 'origin_label') for elt in factors):
            data['origin_label'] = ','.join(elt.origin_label for elt in factors)
        else:
            data['origin_label'] = None

        if all(getattr(elt, 'rational', False) for elt in factors):
            data['rational'] = True
            data['conjugate'] = None

        return cls(data, **kwargs)

    @lazy_attribute
    def origin(self):
        pass

    @lazy_attribute
    def poles(self):
        return []

    @lazy_attribute
    def prelabel(self):
        # this also sets:
        # - nu_*
        # - mu_*
        # and updates:
        # - gamma_factors

        def CCtuple(z):
            return (z.real(), z.imag().abs(), z.imag())

        def spectral_str(x, conjugate=False):
            if conjugate:
                assert x <= 0
                x = -x
                res = "c"
            elif x < 0:
                x = -x
                res = "m"
            else:
                res = "p"
            if x == 0:
                res += "0"
            else:
                res += "%.2f" % x
            return res

        GR, GC = self.gamma_factors
        GR = [CDF(elt) + self.analytic_normalization for elt in GR]
        GC = [CDF(elt) + self.analytic_normalization for elt in GC]
        b, e = Integer(self.conductor).perfect_power()
        if e == 1:
            conductor = b
        else:
            conductor = "{}e{}".format(b, e)
        beginning = "-".join(map(str, [self.degree, conductor, self.central_character]))

        GRcount = Counter(GR)
        GCcount = Counter(GC)
        # convert gamma_R to gamma_C
        for x in sorted(GRcount):
            shift = min(GRcount[x], GRcount[x + 1])
            if shift:
                # We store the real parts of nu doubled
                GCcount[x] += shift
                GRcount[x] -= shift
                GRcount[x] -= shift
        GR = sum([[m] * c for m, c in GRcount.items()], [])
        GC = sum([[m] * c for m, c in GCcount.items()], [])
        assert self.degree == len(GR) + 2 * len(GC)
        GR.sort(key=CCtuple)
        GC.sort(key=CCtuple)

        self.mu_imag = [elt.imag() for elt in GR]
        self.nu_imag = [elt.imag() for elt in GC]

        # deal with real parts
        GR_real = [elt.real() for elt in GR]
        GC_real = [elt.real() for elt in GC]
        self.mu_real = [x.round() for x in GR_real]
        assert set(self.mu_real).issubset({0, 1})
        self.nu_real_doubled = [(2 * x).round() for x in GC_real]
        GRcount = Counter(GR_real)
        GCcount = Counter(GC_real)
        ge = GCD(GCD(list(GRcount.values())), GCD(list(GCcount.values())))
        if ge > 1:
            GR_real = sum(([k] * (v // ge) for k, v in GRcount.items()), [])
            GC_real = sum(([k] * (v // ge) for k, v in GCcount.items()), [])

        rs = ''.join('r%d' % elt.real().round() for elt in GR_real)
        cs = ''.join('c%d' % (elt.real() * 2).round() for elt in GC_real)
        gammas = "-" + rs + cs
        if ge > 1:
            gammas += "e%d" % ge
        if self.algebraic:
            end = "-0"
        else:
            end = "-"
            for G in [GR, GC]:
                for i, elt in enumerate(G):
                    conjugate = False
                    if elt.imag() <= 0 and i < len(G) - 1 and elt.conjugate() == G[i + 1]:
                        conjugate = True
                    elif elt.imag() >= 0 and i > 0 and elt.conjugate() == G[i - 1]:
                        # we already listed this one as a conjugate
                        continue
                    end += spectral_str(elt.imag(), conjugate=conjugate)
        # these are stored in algebraic normalization
        self.gamma_factors = [[elt - self.analytic_normalization for elt in G]
                              for G in [GR, GC]]
        self.spectral_label = gammas + end
        return beginning + gammas + end

    @lazy_attribute
    def st_group(self):
        return None

    @lazy_attribute
    def Lhash(self):
        if self.primitive:
            # this is the integer x st
            # z1 \in [x-1, x+1]*2^-100
            return str(mid_point_100bits(self.positive_zeros_arb[0])[0])
        if hasattr(self, 'factors_obj') and all(getattr(elt, 'Lhash', None) for elt in self.factors_obj):
            return ','.join(elt.Lhash for elt in self.factors_obj)

    @lazy_attribute
    def LhashArray(self):
        if self.primitive:
            return [self.Lhash]
        if hasattr(self, 'factors_obj') and all(getattr(elt, 'Lhash', None) for elt in self.factors_obj):
            return [elt.Lhash for elt in self.factors_obj]

    @lazy_attribute
    def index(self):
        return None

    @lazy_attribute
    def label(self):
        return None

    @lazy_attribute
    def factors(self):
        if self.primitive and hasattr(self, 'label'):
            return [self.label]
        if hasattr(self, 'factors_obj') and all(getattr(elt, 'label', None) for elt in self.factors_obj):
            return sum([elt.factors for elt in self.factors_obj], [])

    @lazy_attribute
    def analytic_normalization(self):
        return 0.5 * self.motivic_weight

    def set_dirichlet(self):
        assert len(self.euler_factors) >= 4
        for i, ai in enumerate(dirichlet_coefficients(self.euler_factors[:4])):
            if i > 1:
                self.__dict__['a' + str(i)] = ai

    # there must be a smarter way to do this
    @lazy_attribute
    def a2(self):
        self.set_dirichlet()
        return self.__dict__['a2']

    @lazy_attribute
    def a3(self):
        self.set_dirichlet()
        return self.__dict__['a3']

    @lazy_attribute
    def a4(self):
        self.set_dirichlet()
        return self.__dict__['a4']

    @lazy_attribute
    def a5(self):
        self.set_dirichlet()
        return self.__dict__['a5']

    @lazy_attribute
    def a6(self):
        self.set_dirichlet()
        return self.__dict__['a6']

    @lazy_attribute
    def a7(self):
        self.set_dirichlet()
        return self.__dict__['a7']

    @lazy_attribute
    def a8(self):
        self.set_dirichlet()
        return self.__dict__['a8']

    @lazy_attribute
    def a9(self):
        self.set_dirichlet()
        return self.__dict__['a9']

    @lazy_attribute
    def a10(self):
        self.set_dirichlet()
        return self.__dict__['a10']

    @lazy_attribute
    def analytic_conductor(self):
        GF_analytic = [tuple(elt + self.analytic_normalization for elt in G)
                       for G in self.gamma_factors]
        return float(self.conductor * conductor_an(*GF_analytic))

    @lazy_attribute
    def root_analytic_conductor(self):
        return float(RDF(self.analytic_conductor).nth_root(self.degree))

    @lazy_attribute
    def bad_primes(self):
        return force_list_int(Integer(self.conductor).prime_divisors())

    @lazy_attribute
    def conductor_radical(self):
        return prod(self.bad_primes)

    @lazy_attribute
    def trace_hash(self):
        if hasattr(self, 'factors_obj') and all(getattr(elt, 'trace_hash', None) for elt in self.factors_obj):
            return sum(elt.trace_hash for elt in self.factors_obj) % 0x1FFFFFFFFFFFFFFF  # mod 2^61 - 1

    @lazy_attribute
    def orbit(self):
        return None

    @lazy_attribute
    def factors_inv(self):
        return []

    @lazy_attribute
    def euler_factors_factorization(self):
        def factorization(original_poly):
            poly = self.ZZT(original_poly)
            assert poly[0] == 1
            if poly == 1:
                return [1]
            facts = poly.factor()
            # if the factor is -1+T^2, replace it by 1-T^2
            # this should happen an even number of times, mod powers
            # if the factor is -1+T^2, replace it by 1-T^2
            # this should happen an even number of times, mod powers
            out = [[-g if g[0] == -1 else g, e] for g, e in facts]
            assert prod(g**e for g, e in out) == poly, "%s != %s" % (prod([g**e] for g, e in out), poly)
            return [[g.list(), e] for g, e in out]
        return force_list_int([factorization(elt) for elt in self.euler_factors])


def tensor_charpoly(f, g):
    R = PolynomialRing(g.parent(), "y")
    y = R.gen()
    A = f(y)
    B = R(g.homogenize(y))
    return B.resultant(A)


def base_change(Lpoly, r):
    R = Lpoly.parent()
    T = R.gen()
    S = PolynomialRing(R, 'u')
    u = S.gen()
    return R(Lpoly(u).resultant(u**r - T))


def sym_pol(L, n):
    Ln = L
    for i in range(2, n + 1):
        Ln = tensor_charpoly(Ln, L)
        b = base_change(L, i)
        extra = (Ln // b).sqrt(2)
        Ln = b * extra
        if Ln[0] == -1:
            Ln *= -1
    return Ln


@cached_function
def sym_pol_gen(n):
    S = PolynomialRing(Integers(), 2, "a, p")
    R = PolynomialRing(S, "T")
    L = R("1 - a*T + p*T^2")
    return sym_pol(L, n).list()


def sym_pol_ECQ(a, p, n):
    return [int(elt.substitute(a=a, p=p)) for elt in sym_pol_gen(n)]


class SmalljacParser(LfunctionsParser):
    ZZT = PolynomialRing(ZZ, "T")
    block_col = 2

    @classmethod
    def ef_multiply_by_riemann(cls, ef, p, power):
        return cls.ZZT(ef) * cls.ZZT([1, -p**(power // 2)])

    def __init__(self, load_key):
        in_headers = in_headers = 'origin_label:power:conductor:curve:hard_factors'.split(':')
        out_headers = 'origin_label:trace_hash:second_moment:root_number_acb:order_of_vanishing:leading_term_arb:special_values_acb:positive_zeros_arb:positive_zeros_extra:plot_delta:plot_values'.split(':')
        self.load_key = load_key
        LfunctionsParser.__init__(self, in_headers, out_headers)

    @lazy_attribute
    def res_constant(self):
        return {
            'coefficient_field': '1.1.1.1',
            'conjugate': None,
            'load_key': self.load_key,
            'self_dual': True,
            'central_character': '1.1',
        }

    def process_record(self, rec):
        if rec.ECcmmod4:
            res = dict(self.res_constant)
            res['euler_factors'] = force_list_int([self.ef_multiply_by_riemann(elt, p, rec.power) for p, elt in zip(primes_first_n(len(rec.euler_factors)), rec.euler_factors)])
            res['bad_lfactors'] = force_list_int([[p, self.ef_multiply_by_riemann(elt, p, rec.power)] for p, elt in rec.bad_lfactors])
            recrh = lfunction_element.from_factors(
                [rec, riemann(rec.power // 2)],
                res
            )
            # move types and urls to the righteous L-func
            recrh.instance_types, recrh.instance_urls = rec.instance_types, rec.instance_urls
            rec.instance_types = rec.instance_urls = []
            recrh.origin = rec.origin
            yield recrh
        yield rec

    def read_curve(self, s):
        if '/' in s:
            # e.g.: [a,0,a,-1,0]/(a^2-23)
            ainv, nfield = s.split('/')
            poly = coeff_to_poly(nfield.strip('()').strip())
            K = NumberField(poly, str(poly.parent().gen()))
            ainv = [K(elt) for elt in ainv.strip('[]').split(',')]
            return EllipticCurve(K, ainv)
        elif 'x' in s:
            # e.g. [-x^2-x, x^3+x^2+x+1]
            s = s.strip()
            if s[0] == '[' and s[-1] == ']':
                R = PolynomialRing(ZZ, 'x')
                f, h = [R(elt) for elt in s.strip('[]').split(',')]
                return HyperellipticCurve(f, h)
            else:
                raise NotImplementedError("only implemented for genus 1 with ainvs")
        else:  # ECQ
            return EllipticCurve(atoii(s))

    def read_hard_factors(self, s):
        primes, factors = s[1:-1].split('],', 1)

        primes += ']'
        primes = atoii(primes, level=1)
        factors = atoii(factors, level=2)
        # FIXME
        # if base field != QQ, some primes might repeated, but in that case we will compute all euler factors directly
        return dict(zip(primes, factors))

    def read_power(self, s):
        return int(s)

    def read_record(self, in_line, out_line):
        res = dict(self.res_constant)
        res.update(self.read_in_line(in_line))
        out_res = self.read_out_line(out_line)
        assert res['origin_label'] == out_res['origin_label'], str((res['origin_label'], out_res['origin_label']))
        res.update(out_res)

        res['motivic_weight'] = res['power']
        res['ECcmmod4'] = False  # overwritten later
        if res['power'] > 1:
            assert res['curve'].genus() == 1
            res['degree'] = res['power'] + 1
            u = ceil(res['power'] / 2)
            res['gamma_factors'] = [([-2 * (u // 2)] if res['power'] % 2 == 0 else []),
                                    sorted([-elt for elt in range(u)])]
            if res['power'] % 4 == 0 and res['curve'].has_cm():
                res['ECcmmod4'] = True
                # account for the missing Riemann zeta factor
                res['degree'] = res['power']
                res['gamma_factors'][0] = []
        else:
            g = res['curve'].genus()
            d = res['curve'].base_ring().degree()
            res['degree'] = 2 * g * d
            res['gamma_factors'] = [[], [0] * g * d]

        # pin down root_number
        assert res['root_number_acb'].contains_integer()
        res['root_number_acb'] = self.CBF(Integer(res['root_number_acb']))
        assert res['root_number_acb'] in [1, -1]
        res['root_angle'] = 0 if res['root_number_acb'] == 1 else 0.5

        return smalljac(res)


class smalljac(lfunction_element):
    def __init__(self, data):
        lfunction_element.__init__(self, data)

    @lazy_attribute
    def euler_factors(self):
        self.set_euler_factors()
        return self.euler_factors

    @lazy_attribute
    def bad_lfactors(self):
        self.set_euler_factors()
        return self.bad_lfactors

    def set_euler_factors(self):
        # Sets:
        # - euler_factors
        # - bad_lfactors
        bound = self.euler_factors_bound
        power = self.power
        bad_primes = Integer(self.conductor).prime_divisors()
        if self.curve.genus() == 1:
            E = self.curve
            T = self.ZZT.gen()
            K = E.base_field()
            if K == QQ:
                N = E.conductor()

                def get_euler_factor(p):
                    Ep = E.local_data(p)
                    if N % p == 0:
                        Lp = 1 - Ep.bad_reduction_type() * T
                    else:
                        Lp = self.ZZT(Ep.minimal_model().change_ring(GF(p)).frobenius_polynomial()).reverse()
                    return Lp

            else:
                N = (E.conductor() * K.discriminant()).absolute_norm()

                def get_euler_factor(p):
                    Lp = 1
                    for f, _ in K.fractional_ideal(p).factor():
                        if f.divides(E.conductor()):
                            local_factor = (1 - E.local_data(f).bad_reduction_type() * T)
                            Lp *= local_factor(T ** f.absolute_norm().valuation(p))
                        else:
                            frob = self.ZZT(E.local_data(f).minimal_model().change_ring(f.residue_field()).frobenius_polynomial())
                            Lp *= frob.reverse()(T ** f.absolute_norm().valuation(p))
                    return Lp
            assert N == self.conductor or power != 1

            if power > 1:
                euler_factors = {}
                for p, v in self.hard_factors.items():
                    euler_factors[p] = v
                for p in prime_range(bound):
                    if p in self.hard_factors:
                        continue
                    ma = get_euler_factor(p)[1]  # -ap
                    euler_factors[p] = sym_pol_ECQ(-ma, p, power)
                if self.ECcmmod4:
                    for p in euler_factors:
                        # euler_factors[p] might be a list
                        euler_factors[p] = self.ZZT(euler_factors[p]) // (1 - T * p**(power // 2))

            else:
                euler_factors = {p: get_euler_factor(p) for p in prime_range(bound) + bad_primes}

        elif self.curve.genus() == 2:
            def get_euler_factor(p):
                if p in self.hard_factors:
                    return self.hard_factors[p]
                else:
                    return self.curve.change_ring(GF(p)).frobenius_polynomial().reverse().list()
            euler_factors = {p: get_euler_factor(p) for p in prime_range(bound) + bad_primes}
        else:
            raise NotImplementedError("only implemented for genus 1 and 2 at the moment")

        self.euler_factors = force_list_int([euler_factors[p] for p in prime_range(bound)])
        self.bad_lfactors = force_list_int([[p, euler_factors[p]] for p in bad_primes])

    @lazy_attribute
    def instance_types(self):
        if self.curve.genus() == 1:
            if self.curve.base_field().degree() == 1:
                if self.power == 1:
                    return ['ECQ']
                else:
                    return ['ECQSymPower']
            else:
                return ['ECNF']
        elif self.curve.genus() == 2:
            return ['G2Q']
        raise NotImplementedError  # another day

    @lazy_attribute
    def origin(self):
        if self.curve.genus() == 1:
            if self.curve.base_field().degree() == 1:
                if self.power == 1:
                    if '^' in self.origin_label:
                        ec_label, power = self.origin_label.split('^')
                        assert int(power) == 1
                    else:
                        ec_label = self.origin_label
                    cond, iso, number = ec_lmfdb_label_regex.match(ec_label).groups()
                    return 'EllipticCurve/Q/%s/%s' % (cond, iso)
                else:
                    ec_label, power = self.origin_label.split('^')
                    assert self.power == int(power)
                    cond, iso, number = ec_lmfdb_label_regex.match(ec_label).groups()
                    return 'SymmetricPower/%d/EllipticCurve/Q/%s/%s' % (self.power, cond, iso)
            else:
                assert len(self.origin_label.split("-")) == 3
                return 'EllipticCurve/' + self.origin_label.replace('-', '/')
        elif self.curve.genus() == 2:
            assert len(self.origin_labe.split(".")) == 2
            return 'Genus2Curve/Q/' + self.origin_label.replace('.', '/')
        raise NotImplementedError  # another day

    @lazy_attribute
    def instance_urls(self):
        return [self.origin]

    @lazy_attribute
    def primitive(self):
        if self.curve.genus() == 1:
            if self.power == 1 and self.curve.base_field().degree() == 1:
                return True
            if self.curve.has_cm() and self.power > 1:
                return False
        return True if self.second_moment == 1 else False

    @lazy_attribute
    def rational(self):
        return True

    @lazy_attribute
    def factors_inv(self):
        if not (self.curve.genus() == 1 and self.curve.has_cm() and self.power > 1):
            return None

        # calls Magma to get the conductor to deduce the various invariants to pindown the factors
        def magma_Lfunction_fac_sympower():
            return magma.Factorisation(magma.SymmetricPower(magma.LSeries(EllipticCurve(list(self.curve.ainvs()))), self.power))

        def magma_Lfunction_inv(L):
            # TODO incorporate trace hash?
            conductor = magma.Conductor(L).sage()
            weight = magma.MotivicWeight(L).sage()
            degree = magma.Degree(L).sage()
            gamma = magma.GammaFactors(L).sage()
            s = str(L)
            translation = 'Translation by '
            if s.startswith(translation):
                shift, _, s = s[len(translation):].split(' ', 2)
                shift = int(shift)
            else:
                shift = 0
            weight = weight - 2 * shift

            gamma = Counter([elt + shift for elt in gamma])
            mus, nus = [], []
            for elt in sorted(gamma):
                if gamma[elt] > 0:
                    if gamma[elt + 1] > 0:
                        gamma[elt] -= 1
                        gamma[elt + 1] -= 1
                        nus.append(2 * elt + weight)
                    else:
                        gamma[elt] -= 1
                        mus.append(elt + weight // 2)

            if 'Grossenchar power' in s:
                for j, elt in [(1, '1st'), (2, '2nd'), (3, '3rd')] + [(i, str(i) + 'th') for i in range(4, 9)]:
                    if s.startswith(elt):
                        assert j == weight
                        break
            if 'Kronecker' in s or 'Riemann zeta function' in s:
                assert weight == 0
                degree = 1
            else:
                degree = 2

            # if we knew the character we could get the prelabel

            return {'label': {'$exists': True},  # we are doing this for the label
                    'degree': degree,
                    'conductor': conductor,
                    'motivic_weight': weight,
                    'rational': True,
                    'mu_real': mus,
                    'mu_imag': [0] * len(mus),
                    'nu_real_doubled': nus,
                    'nu_imag': [0] * len(nus),
                    'rational': True,
                    'primitive': True,
                    'algebraic': True,
                    }

        return [magma_Lfunction_inv(elt[1]) for elt in magma_Lfunction_fac_sympower()
                if 'Riemann zeta function' not in str(elt[1])]


# @cached_function
# def second_moment_to_factors(m):
#    """
#    given the second moment of a representation, returns how it may factor
#    """
#    res = [((1,1),)*m]
#    for k in srange(4, m+1):
#        if k.is_square():
#            for elt in second_moment_to_factors(m-k):
#                res.append(((1,k.sqrt()),) + elt)
#    return res


class artin(lfunction_element):
    def __init__(self):
        raise NotImplementedError


class artin_orbit(artin):
    def __init__(self):
        raise NotImplementedError


class lfunction_collection:
    sep = ':'

    @lazy_class_attribute
    def db(cls):
        from lmfdb import db
        return db

    @lazy_class_attribute
    def tables(cls):
        return [cls.db.lfunc_lfunctions, cls.db.lfunc_data]

    @classmethod
    def dbsearch(cls, query, projection):
        if isinstance(projection, str):
            for table in cls.tables:
                if projection in table.col_type:
                    for elt in table.search(query, projection):
                        yield elt
                else:
                    for _ in range(table.stats._slow_count(query)):
                        yield None
        else:
            assert isinstance(projection, list)
            sp = set(projection)
            for i, table in enumerate(cls.tables):
                for elt in table.search(query,
                                        projection=list(sp.intersection(table.col_type))):
                    yield {k: elt.get(k) for k in sp}

    @lazy_class_attribute
    def lfunctions_schema(cls):
        return OrderedDict((k, cls.db.lfunc_data.col_type[k])
                           for k in sorted(cls.db.lfunc_data.col_type)
                           if k != 'id')

    @lazy_class_attribute
    def instances_schema(cls):
        return OrderedDict((k, cls.db.lfunc_instances.col_type[k])
                           for k in sorted(cls.db.lfunc_instances.col_type)
                           if k != 'id')

    def __init__(self):
        self.lfunctions = defaultdict(list)
        self.ids_to_delete = {}
        self.orbits = defaultdict(list)
        self.instances = []
        self.total = 0

    def populate(self, iterator):
        # this is IO limited
        with Halo(text='Loading collection', spinner='dots') as spinner:
            current = start_time = time.time()
            for i, elt in enumerate(iterator, 1):
                self.lfunctions[elt.prelabel].append(elt)
                if time.time() - current > 1:
                    rate = i / (time.time() - start_time)
                    spinner.text = 'Loading collection: %.2f lines/second' % rate
                    current = time.time()
            old_total = self.total
            self.total = sum(map(len, self.lfunctions.values()))
            spinner.succeed('%d new L-functions loaded in %.2f seconds' % (self.total - old_total, time.time() - start_time))

        def remove_duplicates(objs):
            if len(objs) == 1:
                return objs

            new_objs = objs[:1]
            for x in objs[1:]:
                for y in new_objs:
                    if y == x:
                        y.merge_instances(x)
                        break
                else:
                    # it is a new l-function
                    new_objs.append(x)
            return new_objs

        with Halo(text='Removing duplicates', spinner='dots') as spinner:
            start_time = time.time()
            for prelabel, objs in self.lfunctions.items():
                self.lfunctions[prelabel] = remove_duplicates(objs)
                for elt in self.lfunctions[prelabel]:
                    if elt.orbit:
                        self.orbits[elt.orbit].append(elt)
            old_total = self.total
            self.total = sum(map(len, self.lfunctions.values()))
            spinner.succeed('%d duplicates removed in %.2f seconds' % (old_total - self.total, time.time() - start_time))

    def _populate_label(self):
        info = 'Computing labels'
        with Halo(text=info, spinner='dots') as spinner:
            prelabels_chunks = chunkify([prelabel
                                         for prelabel, objs in self.lfunctions.items()
                                         if any(x.label is None for x in objs)])

            ct = 0
            start_time = time.time()
            for chunk in prelabels_chunks:
                chunk = list(chunk)
                spinner.text = info + progress_bar(ct, self.total, time.time() - start_time)
                db_data = {elt: [] for elt in chunk}
                # spinner.text = info + ': loading data from LMFDB'
                for l in self.dbsearch(
                    {'prelabel': {'$in': chunk}, 'label': {'$exists': True}},
                        projection=lfunction_element.projection + ['origin', 'instance_types', 'instance_urls']):
                    db_data[l['prelabel']].append(l)
                for prelabel, objs in db_data.items():
                    objs = [lfunction_element(l, from_db=True) for l in objs]
                    indexes_taken = [x.index for x in objs] + [x.index for x in self.lfunctions[prelabel] if x.label]
                    nextindex = max(indexes_taken) + 1 if indexes_taken else 0
                    assert nextindex == len(indexes_taken)  # we should not have missing indexes
                    for x in self.lfunctions[prelabel]:
                        ct += 1
                        if x.label:
                            continue
                        for y in objs:
                            if x == y:
                                x.index = y.index
                                x.merge_instances(y)
                                # this distinguishes objects coming from
                                # lfunc_lfunctions vs lfunc_data
                                if getattr(y, 'positive_zeros_mid', None):
                                    self.ids_to_delete[y.id] = y.label
                                break
                        else:
                            x.index = nextindex
                            nextindex += 1

                        x.label = '%s-%s' % (x.prelabel, x.index)

            spinner.succeed('Labels computed in %.2fs' % (time.time() - start_time,))

    def _populate_factors(self):
        projection = [
            'Lhash',
            'accuracy',
            'conductor',
            'degree',
            'label',
            'motivic_weight',
            'order_of_vanishing',
            'positive_zeros',
            'positive_zeros_mid',
            'positive_zeros_rad',
            'primitive',
            'trace_hash',
        ]
        unknown_factors = defaultdict(list)
        unknown_factors_inv = []
        info = 'Computing factors'
        with Halo(text=info, spinner='dots') as spinner:
            for objs in self.lfunctions.values():
                for elt in objs:
                    if elt.primitive or hasattr(elt, 'factors_obj') or elt.factors:
                        continue
                    elif elt.factors_inv:
                        unknown_factors_inv.append(elt)
                    else:
                        unknown_factors[(Integer(elt.conductor), elt.motivic_weight, elt.degree)].append(elt)

            total = len(unknown_factors_inv) + sum(map(len, unknown_factors.values()))
            ct = 0
            start_time = time.time()

            # first go to the ones where we know the invariants
            for obj in unknown_factors_inv:
                ct += 1
                obj.compute_factorization([
                    list(self.dbsearch(query,
                                       projection=projection))
                    for query in obj.factors_inv])

                spinner.text = info + progress_bar(ct, total, time.time() - start_time)

            # this might be quite slow
            for (N, w, d), objs in unknown_factors.items():
                if not N.is_prime():
                    possible_factors = list(self.dbsearch({
                                            'label': {'$exists': True},  # we are doing this for the label
                                            'conductor': {'$in': N.divisors()[1:-1]},
                                            'motivic_weight': {'$lte': w},  # we are also considering translations
                                            'degree': {'$lt': d},
                                            'primitive': True,
                                            'algebraic': True},
                                                          projection=projection))
                    for elt in objs:
                        ct += 1
                        elt.compute_factorization(possible_factors)
                spinner.text = info + progress_bar(ct, total, time.time() - start_time)

            spinner.succeed('Factors computed in %.2fs' % (time.time() - start_time,))

    def orbits_collection(self):
        # generates a new collection from the orbits
        for orbit, objs in self.orbits.items():
            raise NotImplementedError

    def _populate_instances(self):
        info = 'Computing lfun_instances rows'
        with Halo(text=info, spinner='dots') as spinner:
            ct = 0
            start_time = time.time()
            self.instances = []
            lfunctions_dict = {elt.label: elt for objs in self.lfunctions.values() for elt in objs}
            label_chunks = chunkify(lfunctions_dict)
            total = len(label_chunks)
            for chunk in label_chunks:
                ct += 1
                chunk = list(chunk)
                db_data = {elt: set() for elt in chunk}
                for l in self.db.lfunc_instances.search(
                    {'label': {'$in': chunk}},
                        projection=['label', 'type', 'url']):
                    db_data[l['label']].add((l['type'], l['url']))
                for label, pairs in db_data.items():
                    elt = lfunctions_dict[label]
                    orig_pairs = set(zip(elt.instance_types, elt.instance_urls))
                    for t, url in orig_pairs.difference(pairs):
                        self.instances.append({
                            'factors': elt.factors,
                            'label': elt.label,
                            'type': t,
                            'url': url})
                    if pairs:
                        # update elt instance_types and instance_urls
                        elt.instance_types, elt.instance_urls = zip(*orig_pairs.union(pairs))

                spinner.text = info + progress_bar(ct, total, time.time() - start_time)
            spinner.succeed('New lfun_instances rows computed in %.2fs' % (time.time() - start_time,))

    def export(self, lfunctions_filename, instances_filename, overwrite=False, save_memory=False, headers=True):
        if not overwrite and (os.path.exists(lfunctions_filename) or os.path.exists(instances_filename)):
            raise ValueError("Files already exist")

        self._populate_label()
        self._populate_factors()
        self._populate_instances()

        sep = self.sep
        lfunctions_cols = list(self.lfunctions_schema)
        instances_cols = list(self.instances_schema)
        info = 'Writing %s' % lfunctions_filename
        with Halo(text=info, spinner='dots') as spinner:
            ct = 0
            current = start_time = time.time()
            with open(lfunctions_filename, "w") as F:
                if headers:
                    F.write(sep.join(lfunctions_cols))
                    F.write("\n")
                    F.write(sep.join([self.lfunctions_schema[col] for col in lfunctions_cols]))
                    F.write("\n\n")
                ct = 0
                for objs in self.lfunctions.values():
                    ct += len(objs)
                    F.write(
                        '\n'.join(
                            self.sep.join(
                                json_dumps(getattr(elt, col), typ)
                                for col, typ in self.lfunctions_schema.items())
                            for elt in objs) + '\n')
                    if save_memory:
                        objs.clear()
                    if time.time() - current > 1:  # slow down the updates
                        spinner.text = info + progress_bar(ct, self.total, time.time() - start_time)
                        current = time.time()

            spinner.succeed("Wrote %s in %.2f seconds" % (lfunctions_filename, time.time() - start_time))

        info = 'Writing %s' % lfunctions_filename
        with Halo(text=info, spinner='dots') as spinner:
            ct = 0
            current = start_time = time.time()
            total = len(self.instances)
            with open(instances_filename, "w") as F:
                if headers:
                    F.write(sep.join(instances_cols))
                    F.write("\n")
                    F.write(sep.join([self.instances_schema[col] for col in instances_cols]))
                    F.write("\n\n")
                for ct, elt in enumerate(self.instances, 1):
                    # these are dictionaries
                    F.write(sep.join(json_dumps(elt.get(col), self.instances_schema[col]) for col in instances_cols))
                    F.write('\n')
                    if time.time() - current > 1:  # slow down the updates
                        spinner.text = info + progress_bar(ct, total, time.time() - start_time)
                        current = time.time()
            spinner.succeed("Wrote %s in %.2f seconds" % (instances_filename, time.time() - start_time))


class riemann(lfunction_element, CachedRepresentation):
    plot_range = 64
    plot_delta = 1 / 128
    degree = 1
    conductor = 1
    order_of_vanishing = 0
    label = '1-1-1.1-r0-0-0'
    Lhash = '17917892809981029483707922301037'
    positive_zeros_extra = []
    bad_lfactors = []
    self_dual = True
    primitive = True
    rational = True

    def __init__(self, shift):
        # trigger the computation of plot_values
        self.plot_values
        assert(self.positive_zeros_arb[0] in RBF('[14.13472514173469 +/- 5.52e-15]'))
        self.motivic_weight = 2 * shift
        self.gamma_factors = [[-shift], []]
        self.poles = [1 + shift]
        self.origin_label = 'zeta_%s' % shift

    @lazy_class_attribute
    def root_number_acb(cls):
        return cls.CBF(1)

    @lazy_class_attribute
    def leading_term_arb(cls):
        return -cls.RBF.zeta(Integer(1) / 2)

    @lazy_attribute
    def special_values_acb(self):
        return {i: self.RBF.zeta(i) for i in range(2, 10)}

    @lazy_attribute
    def euler_factors(self):
        return [[1, - p**self.shift] for p in prime_range(self.euler_factors_bound)]

    # @lazy_class_attribute
    # def Z(cls):
    #    # documentation |s| = s * sbar
    #    E = -1/8 # (1*(1/2 -1) + 0)/4 # in the documentation is divided by 2
    #    N = 1
    #    gammarabs = lambda s: (gamma(s/2)*s.parent().pi()**(-s/2)).abs()
    #    correction = lambda z: exp(1*z.imag()*z.arg()/2 -E*z.norm().log())*gammarabs(z)*N**(1/4)
    #    return lambda t: cls.zeta.hardy(t)/correction(t.parent().complex_field()(1/2, t))

    @lazy_class_attribute
    def plot_values(cls):
        def theta(t):
            return CDF(0.25, t * 0.5).gamma().arg() - t * CDF.pi().log() / 2

        def Z(t):
            return (CDF(0, theta(t)).exp() * CDF(0.5, t).zeta()).real()

        return [Z(RDF(t)) for t in srange(0, cls.plot_range + cls.plot_delta, cls.plot_delta)]

    @lazy_class_attribute
    def positive_zeros_arb(cls):
        old_prec = gp.set_precision(cls.RBF.precision() * log(2) / log(10) + 10)
        zeros = [cls.RBF(elt) for elt in gp.lfunzeros(gp.lfuninit(1, [13, 79]), [14, 78])]
        gp.set_precision(old_prec)
        return zeros

    @lazy_class_attribute
    def trace_hash(cls):
        from lmfdb.utils.trace_hash import TH_P, TraceHash_from_ap
        return TraceHash_from_ap([1 for _ in TH_P])


def split(fn_in, col, n, sep=LfunctionsParser.sep):
    if n == 1:
        return [(0, -1)]
    num_lines = sum(1 for line in open(fn_in))
    desired_lines = (num_lines + n - 1) // n
    splits = []
    with open(fn_in) as iF:
        searching = False
        oldcol = None
        begin = 0
        for i, line in enumerate(iF):
            if i != 0 and i % desired_lines == 0:
                searching = True
                oldcol = line.split(sep)[col]
            if searching:
                newcol = line.split(sep)[col]
                if newcol != oldcol:
                    splits.append((begin, i))
                    begin = i
                    searching = False
        else:
            splits.append((begin, i + 1))
    return splits


def stitch(main, suffixes):
    if suffixes:
        text = f"Stitching {main} "
        start_time = time.time()
        with Halo(text, spinner="dots") as spinner:
            with open(main, 'a') as W:
                for i, suffix in enumerate(suffixes, 1):
                    filename = main + suffix
                    with open(filename) as F:
                        for line in F:
                            W.write(line)
                    spinner.text = text + progress_bar(i, len(suffixes), time.time() - start_time)
                    os.remove(filename)
            spinner.succeed(f'{main} stitched in in %.2fs' % (time.time() - start_time,))
    return True


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "-j",
        "--jobs",
        dest="jobs",
        metavar="N",
        type=int,
        help="Run up to N jobs in parallel. [default: %(default)s]",
        default=1,
    )
    parser.add_argument(
        "--lmfdb",
        dest="lmfdb",
        metavar="dir",
        help="LMFDB directory [default: %(default)s]",
        default="../lmfdb",
    )
    parser.add_argument(
        "--begin",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--end",
        type=int,
        default=-1
    )
    parser.add_argument('class_name', help='Class name of the reading class, e.g. SmalljacParser')
    parser.add_argument('input', help='file used as input for the lfunctions library')
    parser.add_argument('output', help='the output from the lfunctions library')
    parser.add_argument('lfunc_data', help='file used to add columns to the table lfunc_instances')
    parser.add_argument('lfunc_instances', help='file used to add columns to the table lfunc_instances')
    parser.add_argument('class_args', nargs='*', help='args used to construct the reading class', default=[])

    return parser.parse_args()


def do_block(C, fn_in, fn_out, fn_lfun, fn_ins, begin, end):
    L = lfunction_collection()
    L.populate(C.read_files(fn_in, fn_out, begin, end))
    L.export(fn_lfun, fn_ins, overwrite=True, save_memory=True, headers=(begin == 0))
    return True


def do_block_Popen(args, j, begin, end):
    suffix = '' if begin == 0 else ('_%d' % begin)
    args = ['sage', '-python', os.path.basename(__file__),
            '--lmfdb', args.lmfdb,
            '--begin', str(begin),
            '--end', str(end),
            args.class_name,
            args.input,
            args.output,
            args.lfunc_data + suffix,
            args.lfunc_instances + suffix] + args.class_args

    print(' '.join(args))
    return Halo_wrap_Popen(args, j)


def Halo_wrap_Popen(args, j):
    color = ['grey', 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white'][j % 8]
    prefix = f'{j}:'

    def non_block_read(output):
        ''' even in a thread, a normal read with block until the buffer is full '''
        fd = output.fileno()
        fl = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
        try:
            return output.read()
        except Exception:
            return ''

    def Spinner(message):
        return Halo(prefix + message, spinner="dots", color=color)

    def spinner_text(spinner, output, message):
        last = ''
        for line in output.rstrip('\n').split('\n'):
            if '✔' in line:
                last = ''
                if not spinner:
                    spinner = Spinner(line)
                    spinner.start()
                spinner.succeed(prefix + line)
                spinner = None
            elif line:
                last = line
        else:
            if last:
                if not spinner:
                    spinner = Spinner(last)
                    spinner.start()
                else:
                    spinner.text = prefix + last
        return spinner

    foo = subprocess.Popen(args,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           universal_newlines=True,
                           bufsize=1,
                           )
    message = "Running " + ' '.join(args)
    spinner = Spinner(message)
    spinner.start()
    while foo.poll() is None:
        spinner = spinner_text(spinner,
                               non_block_read(foo.stdout),
                               message)
        time.sleep(0.16)
    else:
        spinner = spinner_text(spinner,
                               non_block_read(foo.stdout),
                               message)
    if foo.returncode != 0:
        print(foo.communicate()[1])
    return foo.returncode


if __name__ == "__main__":
    args = parse_args()
    os.sys.path.append(args.lmfdb)
    stderr = sys.stderr
    f = tempfile.TemporaryFile()
    sys.stderr = open(os.devnull, 'w')
    try:
        from lmfdb import db
        db.logger.propagate = False
        assert db
    except Exception as e:
        sys.stderr = stderr
        raise e
    sys.stderr = stderr
    C = globals()[args.class_name](*args.class_args)

    if args.jobs == 1:
        do_block(C, args.input, args.output,
                 args.lfunc_data, args.lfunc_instances, args.begin, args.end)

    else:
        assert (args.begin, args.end) == (0, -1)
        print('Running %d jobs' % args.jobs)
        splits = split(args.input, C.block_col, args.jobs)
        jobs = []
        with ThreadPoolExecutor(max_workers=args.jobs) as e:
            for j, (begin, end) in enumerate(splits, 1):
                jobs.append(e.submit(do_block_Popen, args, j, begin, end))
        for i, j in enumerate(jobs):
            if j.result() != 0:
                raise ValueError(f'Something went wrong with {i} Exit code:{j.result()}')

        jobs = []
        suffixes = [f'_{elt[0]}' for elt in splits[1:]]
        with ThreadPoolExecutor(max_workers=args.jobs) as e:
            jobs.append(e.submit(stitch, args.lfunc_data, suffixes))
            jobs.append(e.submit(stitch, args.lfunc_instances, suffixes))

        for i, j in enumerate(jobs):
            if not j.result():
                raise ValueError(f'Something went wrong with {i}')
