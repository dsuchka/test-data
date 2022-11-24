#!/usr/bin/python -i

BASE2_LETTERS = '01'
BASE8_LETTERS = '01234567'
BASE10_LETTERS = '0123456789'
BASE16_LETTERS = '0123456789ABCDEF'
BASE32_LETTERS = (
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "234567")
BASE64_LETTERS = (
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/")
BASE58_LETTERS = (
    "  1 2 3 4 5 6 7 8 9"
    "A B C D E F G H   J K L M N   P Q R S T U V W X Y Z"
    "a b c d e f g h i j k   m n o p q r s t u v w x y z").replace(" ", "")

def gcd(a, b):
    '''Returns Greatest Common Divisor of `a` and `b`'''
    while b:
        q = a // b
        (a, b) = (b, a - q * b)
    return abs(a)

def lcm(a, b):
    '''Returns Least Common Multiple of `a` and `b`'''
    return abs(a) * (abs(b) // gcd(a, b))

def extended_gcd(a, b):
    (s, old_s) = (0, 1)
    (t, old_t) = (1, 0)
    (r, old_r) = (b, a)
    while r:
        q = old_r // r
        (old_r, r) = (r, old_r - q * r)
        (old_s, s) = (s, old_s - q * s)
        (old_t, t) = (t, old_t - q * t)
    return (old_r, (old_s, old_t), (t, s))

def inverse_mod(a, p):
    '''Returns `(1/a) mod p`'''
    return extended_gcd(a, p)[1][0] % p

def to_iterable(seq):
    if hasattr(seq, '__next__'): return seq
    return iter(seq)

def continuant_generator(terms, init=(0, 1)):
    terms = to_iterable(terms)
    (v_2, v_1) = init
    while True:
        x = next(terms)
        (v_2, v_1) = (v_1, v_1 * x + v_2)
        yield v_1

def continuant(terms, init=(0, 1), n=None):
    cg = continuant_generator(terms, init=init)
    v = init[1]
    i = -1
    try:
        while (n is None) or (i < n):
            v = next(cg)
            i += 1
    except StopIteration:
        pass
    if (n is not None) and (i < n):
        raise Exception("not enough elements in terms: n={}, terms produced only {}".format(n, i))
    return (v, i)

def const_generator(x):
    while True:
        yield x

def prefixed_generator(prefix_seq, general_seq):
    prefix_seq = to_iterable(prefix_seq)
    try:
        while True:
            yield next(prefix_seq)
    except StopIteration:
        pass
    while True:
        yield next(general_seq)

def prefixed_const_generator(prefix_seq, x):
    return prefixed_generator(prefix_seq, const_generator(x))

def sequential_generator(func, init=1):
    while True:
        yield func(init)
        init += 1

def prefixed_sequential_generator(prefix_seq, func, init=1):
    return prefixed_generator(prefix_seq, sequential_generator(func, init=init))

def repeating_generator(seq):
    seq = to_iterable(seq)
    v = []
    try:
        while True:
            x = next(seq)
            v.append(x)
            yield x
    except StopIteration:
        pass
    if not v: raise StopIteration()
    i = 0
    while True:
        yield v[i]
        i = (i + 1) % len(v)

def prefixed_repeating_generator(prefix_seq, general_seq):
    return prefixed_generator(prefix_seq, repeating_generator(general_seq))

def scf_convergent_generator(terms):
    '''
    Returns a generator that produces successive convergents of a simple continued fraction
      formed by applying the fundamental recurrence formulas:
        * P(n) = P(n-1)*a(n) + P(n-2); P(-2) = 0, P(-1) = 1; [P(0) = a0]
        * Q(n) = Q(n-1)*a(n) + Q(n-2); Q(-2) = 1, Q(-1) = 0; [Q(0) = 1]
      where P(n)/Q(n) is a n-th convergent and equals to the simple continued fraction:
        P(n)/Q(n) = a0 + 1/(a1 + 1/(a2 + 1/(a3 + ...)))

    The `terms` parameter is a sequence (generator or any iterable object) of a(n) values.

    NOTE:
      * `n` starts from 0
      * a0 (first element from the `terms`) is an integer part [i. e. floor(P(n)/Q(n))]
      * `terms` is usually written as: [a0; a1, a2, ...]

    Returned value is a tuple: (P(n), Q(n)).
    '''
    terms = to_iterable(terms)
    (p_2, p_1) = (0, 1)
    (q_2, q_1) = (1, 0)
    while True:
        a = next(terms)
        (p_2, p_1) = (p_1, p_1 * a + p_2)
        (q_2, q_1) = (q_1, q_1 * a + q_2)
        yield (p_1, q_1)

def scf_convergent(terms, n=None):
    '''
    Returns a n-th convergent produced by the `scf_convergent_generator(terms)`.
    :param: x xyu
    * If the `terms` parameter has less that `n` elements, the Exception will be thrown.
    * If the `n` parameter is None and `terms` is an endless sequence, the stack will be
      overflowed and a corresponding exception will be thrown.

    Returned value is a tuple: (P(n), Q(n), n).
    '''
    cg = scf_convergent_generator(terms)
    (p, q) = (1, 0)
    i = -1
    try:
        while (n is None) or (i < n):
            (p, q) = next(cg)
            i += 1
    except StopIteration:
        pass
    if (n is not None) and (i < n):
        raise Exception("not enough elements in terms: n={}, terms produced only {}".format(n, i))
    return (p, q, i)

def gcf_convergent_generator(denom_terms, num_terms):
    '''
    Returns a generator that produces successive convergents of a generalized continued fraction
      formed by applying the fundamental recurrence formulas:
        * P(n) = P(n-1)*a(n) + P(n-2)*b(n); P(-2) = 0, P(-1) = 1; [P(0) = a(0)]
        * Q(n) = Q(n-1)*a(n) + Q(n-2)*b(n); Q(-2) = 1, Q(-1) = 0; [Q(0) = b(0) = 1]
      where P(n)/Q(n) is a n-th convergent and equals to the generalized continued fraction:
        P(n)/Q(n) = a0 + b1/(a1 + b2/(a2 + b3/(a3 + ... + b(n)/a(n))))

    The `demon_terms` parameter is a sequence (generator or any iterable object)
      of a(n) values [n = 0, 1, 2, ...] (so called `partial denominators`).
    The `num_terms` parameter is a sequence (generator or any iterable object)
      of b(n) values [n = 1, 2, 3, ...] (so called `partial numerators`).

    NOTE:
      * `n` starts from 0 for denominators
      * `n` starts from 1 for numerators, so b0 = 1 and
         should not be produced by the `denom_terms`

    Returned value is a tuple: (P(n), Q(n)).
    '''
    denom_terms = to_iterable(denom_terms)
    num_terms = to_iterable(num_terms)
    a = next(denom_terms)
    b = 1
    (p_2, p_1) = (1, a)
    (q_2, q_1) = (0, b)
    yield (p_1, q_1)
    while True:
        a = next(denom_terms)
        b = next(num_terms)
        (p_2, p_1) = (p_1, p_1 * a + p_2 * b)
        (q_2, q_1) = (q_1, q_1 * a + q_2 * b)
        yield (p_1, q_1)

class Quotient:
    def __init__(self, p, q=1):
        if type(p) != int:
            raise ValueError('p is not integer')
        if type(q) != int:
            raise ValueError('q is not integer')
        if q < 1:
            raise ValueError('q is not natural')
        self.p = p
        self.q = q
        self.canonize()
    def canonize(self):
        d = gcd(abs(self.p), self.q)
        if d > 1:
            self.p //= d
            self.q //= d
        return self
    def __repr__(self):
        return "Quotient({}/{})".format(self.p, self.q)
    @staticmethod
    def toQTuple(that):
        if isinstance(that, Quotient):
            return (that.p, that.q)
        elif type(that) == int:
            return (that, 1)
        else:
            raise ValueError('unsupported type: "{}" (value: {})'.format(type(that), that))
    def __iadd0__(self, xp, xq):
        if self.q != xq:
            self.p = (self.p * xq) + (xp * self.q)
            self.q *= xq
        else:
            self.p += xp
        return self.canonize()
    def __iadd__(self, that):
        (xp, xq) = Quotient.toQTuple(that)
        return self.__iadd0__(xp, xq)
    def __add__(self, that):
        return Quotient(self.p, self.q).__iadd__(that)
    def __isub__(self, that):
        (xp, xq) = Quotient.toQTuple(that)
        self.__iadd0__(-xp, xq)
        return self
    def __sub__(self, that):
        return Quotient(self.p, self.q).__isub__(that)
    def __imul0__(self, xp, xq):
        self.p *= xp
        self.q *= xq
        return self.canonize()
    def __imul__(self, that):
        (xp, xq) = Quotient.toQTuple(that)
        return self.__imul0__(xp, xq)
    def __mul__(self, that):
        return Quotient(self.p, self.q).__imul__(that)
    def __itruediv__(self, that):
        (ixq, ixp) = Quotient.toQTuple(that)
        if (ixq == 0): raise ZeroDivisionError("denominator is zero for inverted quotient: {}".format(that))
        if (ixq < 0): (ixp, ixq) = (-ixp, -ixq)
        return self.__imul0__(ixp, ixq) # invert quotien & mul instead of div: x / (a/b) => x * (b/a)
    def __truediv__(self, that):
        return Quotient(self.p, self.q).__itruediv__(that)
    def __abs__(self):
        return Quotient(abs(self.p), self.q)
    def __lt__(self, that):
        if id(self) == id(that): return False
        (xp, xq) = Quotient.toQTuple(that)
        return self.p*xq < xp*self.q
    def __le__(self, that):
        if id(self) == id(that): return True
        (xp, xq) = Quotient.toQTuple(that)
        return self.p*xq <= xp*self.q
    def __eq__(self, that):
        if id(self) == id(that): return True
        (xp, xq) = Quotient.toQTuple(that)
        return (self.p == xp) and (self.q == xq)
    def __ne__(self, that):
        return not self.__eq__(that)
    def __gt__(self, that):
        if id(self) == id(that): return False
        (xp, xq) = Quotient.toQTuple(that)
        return self.p*xq > xp*self.q
    def __ge__(self, that):
        if id(self) == id(that): return True
        (xp, xq) = Quotient.toQTuple(that)
        return self.p*xq >= xp*self.q
    def __pow__(self, n):
        if type(n) != int:
            raise ValueError('n is not integer')
        return Quotient(self.p**n, self.q**n)
    def toPositionalFractionGenerator(self, radix=10):
        (p, q) = (abs(self.p), self.q)
        x = p // q
        p -= q * x

        # 0. produce sign part as None if value is negative
        if self.p < 0:
            yield None

        # 1. produce integer part as: a(m), a(m-1), a(m-2), ..., a(0)
        m = 0
        while x >= (radix**(m+1)):
            m += 1
        while m >= 0:
            y, x = divmod(x, (radix**m))
            m -= 1
            yield y

        # 2. produce radix point as None
        yield None

        # 3. produce fraction part as: b(0), b(1), ...
        while True:
            p *= radix
            x = p // q
            p -= q * x
            yield x
    def toPositionalFraction(self, radixLetters, fractionPrecision=12):
        g = self.toPositionalFractionGenerator(radix=len(radixLetters))
        v = []
        i = next(g)
        if i is None:
            v.append('-')
            i = next(g)
        while i is not None:
            v.append(radixLetters[i])
            i = next(g)
        if fractionPrecision > 0:
            v.append('.')
            while fractionPrecision > 0:
                v.append(radixLetters[next(g)])
                fractionPrecision -= 1
        return ''.join(v)
    def toDecimalFraction(self, fractionPrecision=12):
        return self.toPositionalFraction(BASE10_LETTERS, fractionPrecision=fractionPrecision)
    def toHexadecimalFraction(self, fractionPrecision=12):
        return self.toPositionalFraction(BASE16_LETTERS, fractionPrecision=fractionPrecision)
    def toBase58Fraction(self, fractionPrecision=12):
        return self.toPositionalFraction(BASE58_LETTERS, fractionPrecision=fractionPrecision)
    @staticmethod
    def fromScfConvergent(terms, n=None):
        v = scf_convergent(terms, n)
        return Quotient(v[0], v[1])
    @staticmethod
    def fromScfConvergentWithPrecision(terms, precisionDenominator):
        cg = scf_convergent_generator(terms)
        (p, q) = (1, 0)
        i = -1
        try:
            while precisionDenominator > q*q:
                (p, q) = next(cg)
                i += 1
        except StopIteration:
            pass
        return (Quotient(p, q), i)
