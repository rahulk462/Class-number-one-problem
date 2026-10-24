"""Class numbers of imaginary quadratic fields via reduced binary quadratic forms.

A primitive positive-definite binary quadratic form ax^2 + bxy + cy^2 of
discriminant D = b^2 - 4ac < 0 is *reduced* when

    |b| <= a <= c,  and  b >= 0 whenever |b| == a or a == c.

Every SL2(Z)-equivalence class contains exactly one reduced form, so the class
number h(D) is the number of primitive reduced forms of discriminant D. Since a
reduced form has 3a^2 <= |D|, the search is finite.

For a squarefree d > 0 the ring of integers of Q(sqrt(-d)) has discriminant
D = -d if d = 3 (mod 4) and D = -4d otherwise, and h(Q(sqrt(-d))) = h(D).
"""
from math import gcd, isqrt


def is_squarefree(n: int) -> bool:
    if n < 1:
        return False
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False
        p += 1
    return True


def field_discriminant(d: int) -> int:
    """Discriminant of Q(sqrt(-d)) for squarefree d > 0."""
    if not is_squarefree(d):
        raise ValueError(f"{d} is not squarefree")
    return -d if d % 4 == 3 else -4 * d


def reduced_forms(D: int):
    """List the primitive reduced forms (a, b, c) of discriminant D < 0."""
    if D >= 0 or D % 4 not in (0, 1):
        raise ValueError(f"{D} is not a negative discriminant")
    forms = []
    a_max = isqrt(-D // 3)
    for a in range(1, a_max + 1):
        for b in range(-a + 1, a + 1):  # -a < b <= a
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0):
                continue
            if gcd(gcd(a, b), c) == 1:
                forms.append((a, b, c))
    return forms


def class_number(D: int) -> int:
    """h(D): number of primitive reduced forms of discriminant D < 0."""
    return len(reduced_forms(D))


def kronecker(D: int, n: int) -> int:
    """Kronecker symbol (D/n) for n >= 1."""
    result = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0:
            return 0
        if D % 8 in (3, 5):
            result = -result
    # Jacobi symbol (D/n) for odd n
    a = D % n if n > 1 else 0
    if n == 1:
        return result
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def class_number_analytic(D: int) -> int:
    """Independent check via Dirichlet's class number formula for fundamental D < 0:

        h(D) = -(w / (2|D|)) * sum_{a=1}^{|D|} (D/a) * a,

    with w = 6, 4, 2 roots of unity for D = -3, -4, other.
    """
    w = {-3: 6, -4: 4}.get(D, 2)
    s = sum(kronecker(D, a) * a for a in range(1, -D + 1))
    h, r = divmod(-w * s, 2 * -D)
    assert r == 0
    return h


def class_number_one_fields(limit: int):
    """Squarefree d <= limit with h(Q(sqrt(-d))) = 1."""
    return [d for d in range(1, limit + 1)
            if is_squarefree(d) and class_number(field_discriminant(d)) == 1]


if __name__ == "__main__":
    import sys
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
    print("Class numbers of Q(sqrt(-d)) for small squarefree d:")
    for d in [1, 2, 3, 5, 6, 7, 10, 11, 14, 15, 19, 23, 43, 47, 67, 163]:
        D = field_discriminant(d)
        print(f"  d = {d:4d}  D = {D:5d}  h = {class_number(D)}  forms = {reduced_forms(D)}")
    found = class_number_one_fields(limit)
    print(f"\nSquarefree d <= {limit} with class number 1: {found}")
    print("Nine Heegner numbers recovered:",
          found == [1, 2, 3, 7, 11, 19, 43, 67, 163])
