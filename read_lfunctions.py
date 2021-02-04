import ast
import time
import json
import os
from collections import Counter, defaultdict
from sage.all import (
    CDF,
    ComplexBallField,
    ComplexField,
    EllipticCurve,
    GCD,
    Integer,
    Integers,
    NumberField,
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    RR,
    RealBallField,
    RealIntervalField,
    RealField,
    ZZ,
    cached_function,
    floor,
    gcd,
    lazy_attribute,
    next_prime,
    prime_powers,
    prime_range,
    primes_first_n,
    prod,
    psi,
    vector,
    ceil,
    log,
)
import re
from dirichlet_conrey import DirichletGroup_conrey, DirichletCharacter_conrey


lhash_regex = re.compile(r'^\d+([,]\d+)*$')
ec_lmfdb_label_regex = re.compile(r'(\d+)\.([a-z]+)(\d*)')

# READ utils
def strip(elt):
    return elt.replace(' ', '').strip()

def atoii(elt, level=1, stripped=True):
    s = strip(elt) if not stripped else elt
    assert level >= 1
    if set('[]') == set(s):
        return ast.literal_eval(s)
    assert s[:level] == '['*level and s[-level:] == ']'*level, s
    if level==1:
        return list(map(int, s[1:-1].split(',')))
    else:
        s = s[level:-1].split('['*(level-1))
        return [atoii('['*(level-1) + c.rstrip(','), level=level-1) for c in s]

# WRITE utils
curly_braces = str.translate({'[':'{',']':'}'})
def json_hack(elt, typ):
    if isinstance(elt, str):
        # assuming that we already took care of the right conversion
        return elt
    else:
        res = json.dumps(elt)
        if typ.endswith('[]'):
            res = res.translate(curly_braces)
        return res.replace('null', r'\N')


def realball_to_mid_rad_str(elt, extra_digits=9, max_rad=2**-103):
    # conversion back and forth from binary to decimal may incur various losses
    # 9 extra digits seems to make the trick
    if elt.mid() != 0:
        digits = ceil(log(elt.mid().abs()/(2*elt.rad()))/log(10)) + extra_digits
        mid = elt.mid().str(digits=digits)
    else:
        mid = '0'
    rad = elt.rad().str()
    assert elt in elt.parent()(mid, float(rad))
    return mid, rad


def complexball_to_mid_rad_str(elt, extra_digits=9):
    m1, r1 = realball_to_mid_rad_str(elt.real())
    m2, r2 = realball_to_mid_rad_str(elt.imag())
    return '{%s,%s}' % (m1, m2), '{%s,%s}' % (r1,r2)

# PROPERTIES utils
def extend_multiplicatively(Z):
    for pp in prime_powers(len(Z)-1):
        for k in range(1, (len(Z) - 1)//pp + 1):
            if gcd(k, pp) == 1:
                Z[pp*k] = Z[pp]*Z[k]

def dirichlet_coefficients(euler_factors):
    R = vector(sum(euler_factors), []).base_ring()
    PS = PowerSeriesRing(R)
    pef = list(zip(primes_first_n(len(euler_factors)), euler_factors))
    an_list_bound = next_prime(pef[-1][0])
    res = [1]*an_list_bound
    for p, ef in pef:
        k = RR(an_list_bound).log(p).floor()+1
        foo = (1/PS(ef)).padded_list(k)
        for i in range(1, k):
            res[p**i] = foo[i]
    extend_multiplicatively(res)
    return res



def enough_digits(elt, digits=6):
    return elt.abs().rad_as_ball()*10**digits - 1 < 0


def mid_point_100bits(elt):
    assert elt.rad().log(2) <= -100
    elt100 = elt*2**100
    eltint = elt100.mid().round()
    ball = elt.parent()([eltint -1, eltint + 1])*2**-100
    assert elt in ball
    return eltint, ball

def ball_from_midpoint(elt):
    assert isinstance(elt, str)
    rn = RealNumber(elt)
    eltint = (rn*2**100).round()
    ball = RealBallField(rn.prec())([eltint -1, eltint + 1])*2**-100
    return ball

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
    assert isinstance(elt, str)
    rn = RealNumber(elt)
    ball = RealBallField(rn.prec())(rn)
    assert '.' in elt
    new_rad = 8*10**(-(len(elt.split('.')[-1]) - 1))
    return ball.add_error(new_rad - ball.rad())


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
    arg = arg_hack(foo)/(2*foo.parent().pi())
    while arg > 0.5:
        arg -= 1
    while arg >= -0.5:
        arg += 1
    return arg




@cached_function
def DirGroup(m):
    return DirichletGroup_conrey(m)


@cached_function
def primitivize(label):
    m, n = [ZZ(a) for a in label.split(".")]
    char = DirichletCharacter_conrey(DirGroup(m), n).primitive_character()
    return "%d.%d" % (char.modulus(), char.number())




logpi = RR.pi().log()
log2pi = (2 * RR.pi()).log()

def log_L_inf(s, mu, nu):
    return ((sum([psi((s + elt) / 2) for elt in  mu]) -
             len(mu) * logpi) / 2 +
            sum([psi(s + elt) for elt in nu]) -
            len(nu) * log2pi)

@cached_function
def conductor_an(GR, GC):
    return (2*log_L_inf(1/2, GR, GC).real()).exp()


class read_lfunctions:
    def __init__(self, in_headers, out_headers):
        default_prec = 300
        self.CBF = ComplexBallField(default_prec)
        self.RBF = RealBallField(default_prec)
        self.out_headers = out_headers
        self.in_headers = in_headers

    def from_arb2(self, lower, upper, exp):
        return self.from_acb2(lower, upper, exp, 0, 0, 1)

    def from_acb2(self, lower_real, upper_real, exp_real, lower_imag, upper_imag, exp_imag):
        return self.CBF(self.RBF([lower_real, upper_real])*2**exp_real, self.RBF([lower_imag, upper_imag])*2**exp_imag)

    def read_vector_float(self, elt):
        if elt == '[]':
            return []
        else:
            return [float(x) for x in elt[1:-1].split(',')]
    def read_origin_label(self, elt):
        return elt


    def read_root_number_acb(self, elt):
        return self.from_acb2(*sum(atoii(elt, level=2),[]))

    def read_order_of_vanishing(self, elt):
        return int(elt)

    def read_leading_term_arb(self, elt):
        return self.from_arb2(*atoii(elt))

    def read_special_values_acb(self, elt):
        return [self.from_acb2(*sum(x,[])) for x in atoii(elt, level=3)]

    def read_positive_zeros_arb(self, elt):
        return [self.from_arb2(*x) for x in atoii(elt, level=2)]

    def read_positive_zeros_extra(self, elt):
        return self.read_vector_float(elt)

    def read_plot_delta(self, elt):
        return float(elt)

    def read_plot_values(self, elt):
        return self.read_vector_float(elt)

    def read_mus(self, elt):
        return self.read_vector_float(elt)

    def read_conductor(self, elt):
        return int(elt)


    def read_trace_hash(self, elt):
        return int(elt)

    def read_2nd_moment(self, elt):
        return float(elt)


    def read_line(self, elt, headers, sep):
        return {k: getattr(self, 'read_' + k)(v) for k, v in zip(headers, strip(elt).split(sep))}


    def read_out_line(self, line, sep=":"):
        return self.read_line(line, self.out_headers, sep)

    def read_in_line(self, line, sep=":"):
        return self.read_line(line, self.in_headers, sep)

    def read_record(self, in_line, out_line, sep=":"):
        raise NotImplementedError("you should override this ;)")

    def read_files(self, in_file, out_file):
        with open(in_file) as iF:
            with open(out_file) as oF:
                for i, o in zip(iF, oF):
                    yield self.read_record(i,o)


class lfunction_element:
    projection = [
        'Lhash',
        'conductor',
        'degree',
        'id',
        'index',
        'label',
        'order_of_vanishing',
        'prelabel',
        'primitive',
        'root_angle',
        'trace_hash',
        'positive_zeros',
        'positive_zeros_mid',
        'positive_zeros_rad',
    ]
    def __init__(self, data, from_db=False):
        self.algebraic = True
        self.coeff_info = None
        self.dirichlet_coefficients = None # we prefer ai and euler_factors
        self.credit = None
        self.group = None
        self.st_group = None # we don't have the infrastructure to determine this on the fly
        self.symmetry_type = None
        self.sign_arg = None
        self.from_db = from_db
        if from_db:
            assert set(self.projection).issubset(set(data))
            # convert from literal to integer
            data['conductor'] = int(data['conductor'])

        for i in range(2, 11):
            self.__dict__['A' + str(i)] = None
        self.__dict__.update(data)


        if hasattr(self, 'leading_term_arb'):
            self.leading_term_mid, self.leading_term_rad = realball_to_mid_rad_str(self.leading_term_arb)
            self.leading_term = None

        if hasattr(self, 'root_number_acb'): % 0.5
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

            self.root_number_mid,  self.root_number_rad = complexball_to_mid_rad_str(self.root_number_acb)
            self.sign_arg_mid, self.sign_arg_rad = realball_to_mid_rad_str(arg)
            self.root_angle = float(arg)


        if hasattr(self, 'special_values_acb'):
            if self.dual:
                # pin down the values to the real axis
                for elt in self.special_values_acb:
                    assert elt.imag().contains_zero()
                self.special_values_acb = [elt.parent()(elt.real()) for elt in self.special_values_acb]
            sv = [(i + self.analytic_normalization, complexball_to_mid_rad_str(elt))
                  for (i, elt) in enumerate(self.special_values_acb, 1) if enough_digits(elt)]
            # at the moment artin and smalljac compute at the same points
            self.special_values_at = '{%s}' % ','.join(str(elt[0])for etl in sv)
            self.special_values_mid = '{%s}' % ','.join(elt[1][0] for elt in sv)
            self.special_values_rad = '{%s}' % ','.join(elt[1][1] for elt in sv)
            self.values = None

        if hasattr(self, 'positive_zeros_arb'):
            self.accuracy = self.precision = None # we are using balls
            # we increase the radius to 2**-103 => 31 digits
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

        if from_db:
            # we will use zero balls for comparisons to figure out label and/or factors
            if getattr(self, 'positive_zeros_rad', None) and getattr(self, 'positive_zeros_mid', None):
                R = RealBallField(self.positive_zeros_mid[0].prec())
                self.positive_zeros_arb = [R(m, r) for m, r in zip(self.positive_zeros_mid, self.positive_zeros_rad)]
            elif getattr(self, 'positive_zeros', None):
                if getattr(self, 'accuracy', None) == 100:
                    z1int = self.Lhash.split(',')[0]
                    self.positive_zeros_arb = [ball_from_midpoint(z) for z in self.positive_zeros]
                    assert self.Lhash.split(',')[0] == mid_point_100bits(self.positive_zeros_arb[0])
                else:
                    self.positive_zeros_arb = [numeric_to_ball(z) for z in self.positive_zeros]
            else:
                assert False, '%s' % data





    def __eq__(self, other, pedantic=True):
        # first compare exact invariants
        for attr in ['prelabel', 'primitive', 'order_of_vanishing', 'motivic_weight']:
            if getattr(self, attr) != getattr(other, attr):
                return False



        # non exact invariants

        # check if the first 3 zeros overlap
        for i in range(3):
            if not self.positive_zeros_arb[i].overlaps(positive_zeros_arb[i]):
                return False

        if ((self.root_angle - other.root_angle) % 1).abs() > 1e-7:
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


    def divides(self, other):
        if not (other.conductor % sefl.conductor == 0 and
                self.motivic_weight == other.motivic_weight and
                self.degree <= other.degree and
                self.order_of_vanishing <= other.order_of_vanishing):
            return False

        # check if the first 3 zeros show up in other
        k = 0
        for j, z in enumerate(self.positive_zeros_arb[:3]):
            for i, w in enumerate(other.positive_zeros_arb[k:], k):
                # all elements of w are less than all the elements of z
                if w < z:
                    continue
                # z overlaps with w, we can move to the next zero
                if w.overlaps(z):
                    k = i
                    break
                # all elements of w are greater than all the elements of z
                # and thus it will also be true for every other future w
                if w > z:
                    return False
            else:
                # we got to the end of other.positive_zeros_arb
                # and they were all smaller than z
                assert j > 0 # assuring at least a zero matched
                break
        return True






    def merge_instances(self, other):
        # merge instances
        assert self == other
        assert len(other.instances_types) == len(other.instances_urls)
        self.instances_types += other.instances_types
        self.instances_urls += other.instances_urls
        other.instances_types = self.instances_types
        other.instances_urls = self.instances_urls

    @staticmethod
    def from_factors(factors, data={}):
        """
        one may pass precomputed data via optional data parameter, e.g., euler factors
        """
        assert len(factors) >= 2
        factors.sort(key = lambda elt: (elt.degree, elt.conductor, elt.positive_zeros_arb[0]))
        data['motivic_weight'] = factors[0].motivic_weight
        assert all(data['motivic_weight'] == elt.motivic_weight for elt in factors)
        data['order_of_vanishing'] = sum(elt.order_of_vanishing for elt in factors)
        data['positive_zeros_arb'] = sorted(sum(elt.positive_zeros_arb for elt in factors, []))
        data['conductor'] = prod(elt.conductor for elt in factors)


        if all(hasattr(elt, 'Lhash') for elt in factors):
            data['Lhash'] = ','.join(elt.Lhash)

        if all(hasattr(elt, 'label') for elt in factors):
            data['factors'] = [elt.label for elt in factors]

        if all(hasattr(elt, 'trace_hash') for elt in factors):
            data['trace_hash'] = sum(elt.trace_hash for elt in factors) % 0x1FFFFFFFFFFFFFFF # mod 2^61 -1

        if all(hasattr(elt, 'leading_term_arb') for elt in factors):
            data['leading_term_arb'] = prod(elt.leading_term_arb)

        if all(hasattr(elt, 'root_number_acb') for elt in factors):
            data['root_number_acb'] = prod(elt.leading_term_arbroot_number_acb)
        else:
            # we will use this for comparisons
            arg = sum(elt.root_angle for elt in factors) % 1
            while arg > 0.5:
                arg -= 1
            while arg >= -0.5:
                arg += 1
            data['root_angle'] = arg

        if all(hasattr(elt, 'special_values_acb') for elt in factors):
            data['special_values_acb'] = [prod(sv) for sv in zip(*(elt.special_values_acb for elt in factors))]

        if all(hasattr(elt, 'plot_delta') and hasattr(elt, 'plot_values') for elt in factors):
            factor_plot_values = [ [ ( j * elt.plot_delta,  z) for j, z in enumerate(elt.values) ] for elt in factors]
            interpolations = [spline(elt) for elt in factor_plot_values]
            plot_range = 64.0/data['degree']
            data['plot_delta'] = plot_delta = plot_range/256
            # we cannot hope to get a finer resolution
            assert data['plot_delta'] >= max(elt.plot_delta for elt in factors)
            # we don't want to extrapolate data
            assert all(plot_range <= elt[-1][0] for elt in factor_plot_values)
            data['plot_values'] = [prod([elt(i) for elt in interpolations]) for i in srange(0, plot_range + plot_delta, plot_delta)]
            assert len(data['plot_values']) == 257

        if all(hasattr(elt, 'gamma_factors') for elt in factors):
            data['gamma_factors'] = [sum(elt['gamma_factors'][i] for elt in factors)
                                     for i in range(2)]


        return lfunction_element(data)




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
        GR = [elt + self.analytic_normalization for elt in GR]
        GC = [elt + self.analytic_normalization for elt in GC]
        b, e = self.conductor.perfect_power()
        if e == 1:
            conductor = b
        else:
            conductor = "{}e{}".format(b, e)
        beginning = "-".join(map(str, [self.degree, conductor, self.central_character]))

        GRcount = Counter(GR)
        GCcount = Counter(GC)
        # convert gamma_R to gamma_C
        for x in sorted(GRcount):
            shift = min(GRcount[x], GRcount[x+1])
            if shift:
                # We store the real parts of nu doubled
                GCcount[x] += shift
                GRcount[x] -= shift
                GRcount[x] -= shift
        GR = sum([[m]*c for m, c in GRcount.items()], [])
        GC = sum([[m]*c for m, c in GCcount.items()], [])
        assert self.degree == len(GR) + 2*len(GC)
        GR.sort(key=CCtuple)
        GC.sort(key=CCtuple)

        self.mu_imag = [elt.imag() for elt in GR]
        self.nu_imag = [elt.imag() for elt in GC]

        # deal with real parts
        GR_real = [elt.real() for elt in GR]
        GC_real = [elt.real() for elt in GC]
        self.mu_real = [x.round() for x in GR_real]
        assert set(self.mu_real).issubset(set([0,1]))
        self.nu_real_doubled = [(2*x).round() for x in GC_real]
        GRcount = Counter(GR_real)
        GCcount = Counter(GC_real)
        ge = GCD(GCD(list(GRcount.values())), GCD(list(GCcount.values())))
        if ge > 1:
            GR_real = sum(([k]*(v//ge) for k, v in GRcount.items()), [])
            GC_real = sum(([k]*(v//ge) for k, v in GCcount.items()), [])

        rs = ''.join(['r%d' % elt.real().round() for elt in GR_real])
        cs = ''.join(['c%d' % (elt.real()*2).round() for elt in GC_real])
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
        self.gamma_factors = [[elt -  self.analytic_normalization for elt in G]
                              for G in [GR, GC]]
        self.spectral_label = gammas + end
        return beginning + gammas + end


    @lazy_attribute
    def st_group(self):
        return None

    @lazy_attribute
    def set_primitive(self):
        return None

    @lazy_attribute
    def Lhash(self):
        if self.primitive:
            # this is the integer x st
            # z1 \in [x-1, x+1]*2^-100
            return str(mid_point_100bits(self.zeros_arb[0]))

    @lazy_attribute
    def LhashArray(self):
        if self.primitive:
            return self.Lhash.split(',')

    @lazy_attribute
    def label(self):
        return None

    @lazy_attribute
    def index(self):
        return None

    @lazy_attribute
    def factors(self):
        if self.primitive and hasattr(self, 'label'):
            return [self.label]




    @lazy_attribute
    def analytic_normalization(self):
        return 0.5*self.motivic_weight

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
        GF_analytic = [tuple(elt + self.analytic_normalization) for elt in G)
                       for G in self.gamma_factors]
        return self.conductor * conductor_an(*GF_analytic)


    @lazy_attribute
    def badprimes(self):
        return self.conductor.prime_divisors()


    @lazy_attribute
    def conductor_radical(self):
        return prod(self.bad_primes)


    @lazy_attribute
    def trace_hash(self):
        return None


    @lazy_attribute
    def orbit(self):
        return None
    #@lazy_attribute
    #def instances_urls(self):
    #    self.instances_types = []
    #    return []

    #@lazy_attribute
    #def intances_types(self):
    #    self.instances_urls = []
    #    return []






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
        extra = (Ln//b).sqrt(2)
        Ln = b*extra
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


class read_smalljac(read_lfunctions):
    def __init__(self, load_key):
        in_headers = in_headers ='origin_label:power:conductor:curve:hard_factors'.split(':')
        out_headers = 'origin_label:trace_hash:second_moment:root_number_acb:order_of_vanishing:leading_term_arb:special_values_acb:positive_zeros_arb:positive_zeros_extra:plot_delta:plot_values'.split(':')
        read_lfunctions.__init__(self, in_headers, out_headers)

    @lazy_attribute
    def res_constant(self):
        return {
            'coefficient_field': '1.1.1.1',
            'conjugate': None,
            'load_key': self.load_key,
            'self_dual': True,
            'central_character': '1.1',
        }

    def read_curve(self, s):
        if '/' in s:
            raise NotImplementedError("only implemented for Q")
        elif 'x' in s:
            raise NotImplementedError("only implemented for genus 1 with ainvs")
        else: #ECQ
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

    def read_record(self, in_line, out_line, sep=":"):
        res = dict(self.res_constant)
        res.update(self.read_in_line(in_line))
        out_res = self.read_out_line(out_line)
        assert res['origin_label'] == out_res['origin_label'], str((res['origin_label'], out_res['origin_label']))
        res.update(out_res)

        res['motivic_weight'] = res['power']
        if res['power'] > 1:
            res['degree'] = res['power'] + 1
            u = ceil(res['power']/2)
            res['gamma_factors'] = [([-2*floor(u/2)] if res['power']%2 == 0 else []),
                                    sorted([-elt for elt in range(0, u)])]
        else:
            g = res['curve'].genus()
            res['degree'] = 2*g
            res['gamma_factors'] = [[], [0]*g]

        # pin down root_number
        assert res['root_number_acb'].contains_integer()
        res['root_number_acb'] = self.CBF(Integer(res['root_number_acb']))
        assert res['root_number_acb']  in [1, -1]
        res['root_angle'] = 0 if res['root_number_acb'] == 1 else 0.5


        return smalljac(res)



class smalljac(lfunction_element):
    def __init__(self, data):
        lfunction_element.__init__(self, data)
        self.set_euler_factors()




    def set_euler_factors(self):
        # Sets:
        # - euler_factors
        # - bad_lfactors
        bound = 100
        power = self.power
        if self.curve.genus() == 1:
            E = self.curve
            ZZT = PolynomialRing(ZZ, "T")
            T = ZZT.gen()
            K = E.base_field()
            if K == QQ:
                K = NumberField(T,"a") # making sure SAGE sees K as a numberfield
                E = E.change_ring(K)
            N = (E.conductor() * K.discriminant()**2).absolute_norm()
            assert N == self.conductor


            def get_eulerfactor(p):
                Lp = 1
                for f, _ in K.fractional_ideal(p).factor():
                    f = K.fractional_ideal(f)
                    if f.divides(N):
                        local_factor = (1 - E.local_data(f).bad_reduction_type() * T)
                        Lp *=  local_factor( T ** f.absolute_norm().valuation(p) )
                    else:
                        frob = ZZT(E.local_data(f).minimal_model().change_ring(f.residue_field()).frobenius_polynomial())
                        Lp *= frob.reverse()( T ** f.absolute_norm().valuation(p) )
                return list(map(int, Lp))

            if power > 1:
                euler_factors = {}
                for p, v in self.hard_factors.items():
                    euler_factors[p] = v
                for p in prime_range(bound):
                    if p in self.hard_factors:
                        continue
                    _, ma, p = euler_factors(p)
                    euler_factors[p] = sym_pol_ECQ(-ma, p, power)
            else:
                euler_factors = {p: get_eulerfactor(p) for p in prime_range(bound)}

        else:
            raise NotImplementedError("only implemented for genus 1 at the moment")

        self.euler_factors = [euler_factors[p] for p in prime_range(bound)]
        self.bad_lfactors = [[p, euler_factors[p]] for p in Integer(self.conductor).prime_divisors()]

    @staticmethod
    def euler_factors_factorization(self):
        if self.degree <= 10:
            return None
        raise NotImplementedError # another day

    @lazy_attribute
    def types(self):
        if self.power > 1 and self.curve.genus() == 1:
            return ['SymPower']
        raise NotImplementedError # another day

    @lazy_attribute
    def origin(self):
        if self.curve.genus() == 1 and self.power > 1:
            ec_label, power = self.origin_label.split('^')
            assert self.power == int(power)
            cond, iso, number = ec_lmfdb_label_regex.match(ec_label).groups()
            return 'SymmetricPower/%d/EllipticCurve/Q/%d/%s' % (self.power, cond, iso)
        raise NotImplementedError # another day

    @lazy_attribute
    def instance_types(self):
        if self.curve.genus() == 1 and self.power > 1:
            self.instances_urls = [self.origin]
            return ['SymPower']
        raise NotImplementedError # another day

    @lazy_attribute
    def instances_urls(self):
        if self.curve.genus() == 1 and self.power > 1:
            self.instance_types = ['SymPower']
            return [self.origin]

        raise NotImplementedError # another day

    @lazy_attribute
    def primitive(self):
        return True if round(self.second_moment) == 1 else False


    @lazy_attribute
    def rational(self):
        return True

#@cached_function
#def second_moment_to_factors(m):
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

class lfunction_collection(read_lfunctions):
    def __init__(self, parent, sep=':'):
        from lmfdb import db
        self.db = db
        self.db.lfunc_lfunctions._include_nones
        self.lfunctions_schema = dict(db.lfunc_lfunctions.col_type())
        self.instances_schema = dict(db.lfunc_instances.col_type())
        self.lfunctions_schema.pop('id')
        self.instances_schema.pop('id')
        self.parent = parent
        self.lfunctions = defaultdict(list)
        self.sep = sep
        self.ids_to_delete = {}
        self.orbits = defaultdict(list)





    def populate(self, iterator):
        for elt in iterator:
            self.lfunctions[elt.prelabel].append(elt)

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

        for prelabel, objs in self.lfunctions.items():
            self.lfunctions[prelabel] = remove_duplicates(objs)
            for elt in self.lfunctions[prelabel]:
                if elt.orbit:
                    self.orbits[elt.orbit].append(elt)


    @staticmethod
    def chunkify(l, chunk_size=1000):
        res = [list(l)[elt*size:(elt + 1)*size] for elt in range(len(l)//size + 1)]
        assert len(l) == sum(map(len, res))
        return res


    def _populate_label(self):
        prelabels_chunks = self.chunkify([prelabel for prelabel, objs self.lfunctions.items()
                                          if any(x.label is None in objs)])
        for chunk in prelabels_chunks:
            db_data = defaultdict(list)
            for l in self.db.lfunc_lfunctions.search(
                    {'prelabel': {'$in': chunk}, 'label': {'$exists': True}},
                projection=projection.lfunction_element):
                db_data[l['prelabel']].append(lfunction_element(l, from_db=True))
            for prelabel, objs in db_data.items():
                indexes_taken = [x.index for x in objs] + [x.index for x in self.lfunctions[prelabel] if x.label]
                nextindex = max(indexes_taken) if indexes_taken else 1
                assert nextindex == len(indexes_taken) + 1 # we should not have missing indexes
                for x in self.lfunctions[prelabel]:
                    if x.label:
                        continue
                    for y in objs:
                        if x == y:
                            x.index = y.index
                            x.merge(y)
                            self.ids_to_delete[y.id] = y.label
                            break
                    else:
                        x.index = nextindex
                        nextindex += 1

                    x.label = '%s-%s' % (x.prelabel, x.index)


    def _populate_factors_and_LhashArray(self):
        unknown_factors = defaultdict(list)
        for prelabel, objs in self.lfunctions.items():
            for elt in objs:
                if elt.factors and elt.LhashArray:
                    continue
                else:
                    unknown_factors[(elt.conductor, elt.motivic_weight, elt.degree)].append(elt)

        def write_factors(obj, possible_factors):
            # we are assuming that it factors with no repeated factors
            res = []
            for elt in possible_factors:
                if elt.divides(obj):
                    res.append(obj)

            def error_msg():
                return '%s, %s' % (obj.__dict__, [elt.__dict__ for elt in res])

            if res:
                Lprod = lfunction_element.from_factors(res)
                if obj.degree == Lprod.degree:
                    assert obj.__eq__(Lprod, pedantic=False), error_msg()
                    obj.factors = Lprod.factors
                    obj.Lhash = Lprod.Lhash
                    obj.LhashArray = Lprod.LhashArray
                else:
                    assert Lprod.divides(obj), error_msg()
                    _ = Lprod.prelabel # trigger the computation of the prelabel, so that gamma_factors are in reduced shape
                    for i in range(2):
                        Gprod = Counter(Lprod.gamma_factors[i])
                        Gobj = Counter(obj.gamma_factors[i])
                        for elt, c in Gprod.items():
                            assert Gobj[elt] >= v, error_msg()
                    obj.factors = Lprod.factors + [None]
                    obj.Lhash = "_" + Lprod.Lhash
                    obj.LhashArray = Lprod.LhashArray + [None]
            else:
                obj.factors = [None]
                obj.Lhash = "_" + str(mid_point_100bits(self.zeros_arb[0]))
                obj.LhashArray = [obj.Lhash]






        # this might be quite slow
        for (N, w, d), objs in unknown_factors.items():
            possible_factors = [ lfunction_element(elt, from_db=True) for elt in
                                self.db.lfunc_lfunctions.search(
                                    {
                                        'conductor': {'$in': Integer(N).divisors()},
                                        'motivic_weight': w,
                                        'degree': {'$lt': d},
                                        'primitive': True,
                                        'algebraic': True,
                                    },
                                    projection=lfunction_element.projection + ['gamma_factors'])]
            for elt in obs:
                write_factors(elt, possible_factors)





    def orbits_collection(self):
        # generates a new collection from the orbits
        for orbit, objs in self.orbits.items():
            raise NotImplementedError

    def _instances_rows(self):
        res = []

        for elt in self.lfunctions:
            for t, url in zip(elt.instance_types, elt.instances_urls):
                res.append({
                    'Lhash': elt.Lhash,
                    'LhashArray': elt.LhashArray,
                    'factors': elt.factors,
                    'label': elt.label,
                    'type': t,
                    'url': url})

        return res




    def export(self, lfunctions_filename, instances_filename, overwrite=False):
        if not overwrite and (os.path.exists(lfunctions_filename) or os.path.exists(instances_filename)):
            raise ValueError("Files already exist")


        self._populate_label()
        self._populate_factors_and_LhashArray()
        instances = self._instance_rows()

        print("Writing to %s and %s" % (lfunctions_filename, instances_filename))


        sep = self.sep
        lfunctions_cols = list(self.lfunctions_schema)
        instances_cols = list(self.instances_schema)
        with open(lfunctions_filename, "w") as F:
            start_time = time.time()
            F.write(sep.join(lfunctions_cols))
            F.write("\n")
            F.write(sep.join([self.lfunctions_schema[col] for col in lfunctions_cols]))
            F.write("\n")
            for elt in self.lfunctions:
                F.write(sep.join(json_hack(elt.getattr(col)) for col in lfunctions_cols))
            print("Wrote %s in %.2f seconds" % (lfunctions_filename, time.time() - start_time))

        with open(instances_filename, "w") as F:
            F.write(sep.join(instances_cols))
            F.write("\n")
            F.write(sep.join([self.instances_schema[col] for col in instances_cols]))
            F.write("\n")
            for elt in instances:
                F.write(sep.join(json_hack(elt.getattr(col)) for col in instances_cols))
            print("Wrote %s in %.2f seconds" % (instances_filename, time.time() - start_time))





# TODO


if False:
    desired_schema = {
        'A10': 'numeric',
        'A2': 'numeric',
        'A3': 'numeric',
        'A4': 'numeric',
        'A5': 'numeric',
        'A6': 'numeric',
        'A7': 'numeric',
        'A8': 'numeric',
        'A9': 'numeric',
        'Lhash': 'text',
        'a10': 'jsonb',
        'a2': 'jsonb',
        'a3': 'jsonb',
        'a4': 'jsonb',
        'a5': 'jsonb',
        'a6': 'jsonb',
        'a7': 'jsonb',
        'a8': 'jsonb',
        'a9': 'jsonb',
        'accuracy': 'smallint', # bits of precision for the zeros used by ECNF and CMF
        'algebraic': 'boolean', # = arithmetic
        'analytic_conductor': 'double precision',
        'analytic_normalization': 'numeric', # this could be converted to a double, as it is always an half integer
        'bad_lfactors': 'jsonb',
        'bad_primes': 'bigint[]',
        'central_character': 'text',
        'coeff_info': 'jsonb', # only used by Dirichlet L-functions, as ai and dirichlet_coefficients are stored algebraically
        'coefficient_field': 'text',
        'conductor': 'numeric',
        'conductor_radical': 'integer',
        'conjugate': 'text',
        'credit': 'text',
        'degree': 'smallint',
        'dirichlet_coefficients': 'jsonb',
        'euler_factors': 'jsonb',
        'euler_factors_factorization': 'jsonb', # we only set these for large degree and weight
        'factors': 'text[]', # array with the labels of its factors
        'gamma_factors': 'jsonb', # these are stored in algebraic normalization
        'group': 'text', # \in ['GL2', 'GL3', 'GL4', 'GSp4', None]
        'id': 'bigint', # automatically set
        'index': 'smallint', # set in a second iteration
        'instance_types': 'text[]',
        'instance_urls': 'text[]',
        'label': 'text', # set in a second iteration
        'leading_term': 'text',
        'leading_term_mid': 'numeric',
        'leaving_term_rad': 'real',
        'load_key': 'text',
        'motivic_weight': 'smallint',
        'mu_imag': 'numeric[]',
        'mu_real': 'smallint[]',
        'nu_imag': 'numeric[]',
        'nu_real_doubled': 'smallint[]',
        'order_of_vanishing': 'smallint',
        'origin': 'text',
        'plot_delta': 'numeric',
        'plot_values': 'jsonb',
        'positive_zeros': 'jsonb',
        'positive_zeros_mid': 'numeric[]',
        'positive_zeros_rad': 'double precision[]',
        'positive_zeros_extra': 'double precision[]',
        'precision': 'smallint', # used by MaassForms, it is stored in digits
        'prelabel': 'text',
        'primitive': 'boolean',
        'rational': 'boolean',
        'root_analytic_conductor': 'double precision',
        'root_angle': 'double precision', # stored between -.5 to .5
        'root_number': 'text',
        'root_number_mid': 'numeric[]', # mid real and mid imag
        'root_number_rad': 'numeric[]', # respective radius
        'self_dual': 'boolean',
        'sign_arg': 'numeric',
        'sign_arg_mid': 'numeric',
        'sign_arg_rad': 'real',
        'special_values_at': 'double precision[]', # in algebraic normalization
        'special_values_mid': 'numeric[]', # an array of pairs of mid real and mid imag values
        'special_values_rad': 'real[]', # the respective radius
        'spectral_label': 'text', # copy data
        'st_group': 'text', # we need to set it later
        'symmetry_type': 'text', # \in ['orthogonal', 'symplectic', 'unitary', None]
        'trace_hash': 'bigint',
        'types': 'jsonb',
        'values': 'jsonb', # only used by Dirichlet L-functions
        'z1': 'numeric',
        'z2': 'numeric',
        'z3': 'numeric'
    }





