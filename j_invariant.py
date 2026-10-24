"""The j-invariant at CM points and why e^(pi*sqrt(163)) is nearly an integer.

j(tau) = E4(tau)^3 / Delta(tau) with q = e^(2 pi i tau):

    E4    = 1 + 240 * sum sigma_3(n) q^n
    Delta = q * prod (1 - q^n)^24
    j     = 1/q + 744 + 196884 q + 21493760 q^2 + ...

For tau = (1 + sqrt(-d))/2 we get q = -e^(-pi sqrt d), so

    j(tau) = -e^(pi sqrt d) + 744 - 196884 e^(-pi sqrt d) + ...

When Q(sqrt(-d)) has class number 1, CM theory says j(tau) is a rational
integer. For d = 163 it equals -640320^3, so

    e^(pi sqrt 163) = 640320^3 + 744 - 196884 e^(-pi sqrt 163) + ...

and the error term is about 7.5e-13.
"""
from mpmath import mp, mpf, mpc, sqrt, exp, pi, nint, cbrt, fabs

from class_number import field_discriminant

mp.dps = 60


def j(tau, terms=60):
    q = exp(2j * pi * tau)
    e4 = 1 + 240 * sum(sum(k**3 for k in range(1, n + 1) if n % k == 0) * q**n
                       for n in range(1, terms))
    prod = mpc(1)
    for n in range(1, terms):
        prod *= (1 - q**n)
    delta = q * prod**24
    return e4**3 / delta


def cm_point(d):
    """Generator tau of the ring of integers of Q(sqrt(-d)), in the upper half plane."""
    D = field_discriminant(d)
    return (1 + sqrt(mpf(-D)) * 1j) / 2 if D % 4 == 1 else sqrt(mpf(d)) * 1j


HEEGNER = [1, 2, 3, 7, 11, 19, 43, 67, 163]


if __name__ == "__main__":
    print("j at the CM point of each class-number-one field:")
    for d in HEEGNER:
        val = j(cm_point(d))
        n = int(nint(val.real))
        c = int(nint(cbrt(fabs(n)))) * (1 if n >= 0 else -1)
        print(f"  d = {d:3d}  j = {n:>22d} = ({c})^3   "
              f"|j - round(j)| = {mp.nstr(fabs(val - n), 3)}")

    print("\nFor comparison, d = 5 (class number 2): j =",
          mp.nstr(j(cm_point(5)).real, 25), "(not an integer)")

    x = exp(pi * sqrt(163))
    N = 640320**3 + 744
    print(f"\ne^(pi*sqrt(163))   = {mp.nstr(x, 45)}")
    print(f"640320^3 + 744     = {N}")
    print(f"difference         = {mp.nstr(x - N, 15)}")
    print(f"-196884*e^(-pi*sqrt163) = {mp.nstr(-196884 * exp(-pi * sqrt(163)), 15)}")
    print("within 1e-12 of an integer:", fabs(x - N) < mpf("1e-12"))
