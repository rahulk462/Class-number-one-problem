import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mpmath import mpf, exp, pi, sqrt, fabs, nint

from class_number import (class_number, class_number_analytic,
                          class_number_one_fields, field_discriminant, is_squarefree)
from j_invariant import HEEGNER, cm_point, j


def test_known_class_numbers():
    # h(D) for the fields Q(sqrt(-d)), from standard tables.
    known = {1: 1, 2: 1, 3: 1, 5: 2, 6: 2, 7: 1, 10: 2, 14: 4, 15: 2, 17: 4,
             21: 4, 23: 3, 26: 6, 30: 4, 31: 3, 47: 5, 71: 7, 199: 9, 427: 2}
    for d, h in known.items():
        assert class_number(field_discriminant(d)) == h, d


def test_forms_agree_with_dirichlet_formula():
    for d in range(1, 600):
        if is_squarefree(d):
            D = field_discriminant(d)
            assert class_number(D) == class_number_analytic(D), d


def test_heegner_numbers():
    assert class_number_one_fields(5000) == HEEGNER


def test_j_is_a_cube_integer_at_heegner_points():
    expected = {1: 12**3, 2: 20**3, 3: 0, 7: -15**3, 11: -32**3, 19: -96**3,
                43: -960**3, 67: -5280**3, 163: -640320**3}
    for d, value in expected.items():
        z = j(cm_point(d))
        assert fabs(z - value) < mpf("1e-30"), d


def test_j_not_integer_for_class_number_two():
    z = j(cm_point(5)).real
    assert fabs(z - nint(z)) > mpf("0.01")


def test_ramanujan_constant():
    err = exp(pi * sqrt(163)) - (640320**3 + 744)
    assert mpf("-1e-12") < err < 0
    # The error is the next q-expansion term, -196884 e^(-pi sqrt 163), to ~1e-28.
    assert fabs(err + 196884 * exp(-pi * sqrt(163))) < mpf("1e-25")
