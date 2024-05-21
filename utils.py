import ast
import json
import re
from collections.abc import Iterable
from dirichlet_conrey import DirichletGroup_conrey, DirichletCharacter_conrey
from sage.all import (
    ComplexBallField,
    Integer,
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    RBF,
    RDF,
    RR,
    RealBallField,
    ZZ,
    cached_function,
    gcd,
    next_prime,
    prime_powers,
    primes_first_n,
    prod,
    vector,
)
from sage.rings.real_arb import RealBall
from sage.rings.real_double import RealDoubleElement
from sage.rings.real_mpfr import RealLiteral, RealField, RealNumber


lhash_regex = re.compile(r'^\d+([,]\d+)*$')
ec_lmfdb_label_regex = re.compile(r'(\d+)\.([a-z]+)(\d*)')

def mod1(elt):
    return RDF(elt) - RDF(elt).floor()

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
    as one usually doesn't truncate when printing floats to be able to recover
    all the binary digits

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


def complex_string_to_ball(elt):
    """
    Convert a string representing a complex number into a complex ball by
    setting the radius such that the last two digits of the real parts and
    imaginary parts are unknown.
    """
    assert isinstance(elt, str)
    assert "." in elt   # don't handle trivial strings
    split_symbol = None
    if '+' in elt:
        split_symbol = '+'
    elif '-' in elt:
        split_symbol = '-'
    if split_symbol:
        real_part, imag_part = elt.split(split_symbol)
        sigfig_mantissa_real = real_part.lstrip('-+0.')
        sigfigs_real = len(sigfig_mantissa_real) - ('.' in sigfig_mantissa_real)
        bits_real = int(LOG_TEN_TWO_PLUS_EPSILON * sigfigs_real) + 1
        rad_len_real = len(real_part.split('.')[-1])
        rad_real = 10**(-(rad_len_real - 3))   # assume last 3 digits fuzzy

        sigfig_mantissa_imag = imag_part.lstrip('-+0.').rstrip('iI*')
        sigfigs_imag = len(sigfig_mantissa_imag) - ('.' in sigfig_mantissa_imag)
        bits_imag = int(LOG_TEN_TWO_PLUS_EPSILON * sigfigs_imag) + 1
        rad_len_imag = len(imag_part.split('.')[-1])
        rad_imag = 10**(-(rad_len_imag - 3))   # assume last 3 digits fuzzy

        bits = max(bits_real, bits_imag)
        rad = max(rad_real, rad_imag)

        real_part_ball = RealBallField(bits)(real_part, rad)
        imag_part_ball = RealBallField(bits)(imag_part.rstrip('iI*'), rad)

        return ComplexBallField(bits)(real_part_ball, imag_part_ball)
    else:  # purely real or purely imaginary
        sigfig_mantissa = elt.lstrip('-+0.').rstrip('iI*')
        sigfigs = len(sigfig_mantissa) - ('.' in sigfig_mantissa)
        bits = int(LOG_TEN_TWO_PLUS_EPSILON * sigfigs) + 1
        rad_len = len(elt.rstrip('iI*').split('.')[-1])
        rad = 10**(-(rad_len - 3))   # assume last 3 digits fuzzy
        part_ball = RealBallField(bits)(elt.rstrip('iI*'), rad)
        if "i" in elt or "I" in elt:
            return ComplexBallField(bits)(0, part_ball)
        else:
            return ComplexBallField(bits)(part_ball, 0)


def realnumber_to_ball(elt, R):
    return R(elt, float(elt.ulp()))


def approx_ball(elt, prec=53):
    """
    if we can approximate the ball, returns such approximation
    """
    # this is what we would get from such approximation
    approx_ball = realnumber_to_ball(elt.numerical_approx(prec=prec), RealBallField(prec))
    if elt in approx_ball: # this checks that the ball given is inside of the implicit ball
        return approx_ball.mid()
    else:
        raise RuntimeError("could not approximate ball to the desired precision")


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


