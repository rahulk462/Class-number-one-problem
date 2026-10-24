# Class-number-one-problem

Compute class numbers of imaginary quadratic fields Q(√−d), recover the nine
Heegner numbers (1, 2, 3, 7, 11, 19, 43, 67, 163), and explain why
e^(π√163) is within 10⁻¹² of an integer using the j-invariant.

```
pip install -r requirements.txt
python3 class_number.py [limit]   # class numbers + Heegner search (default limit 10000)
python3 j_invariant.py            # j at CM points and e^(π√163)
python3 -m pytest -q              # checks
```

## 1. Class numbers from binary quadratic forms

For squarefree d > 0 the field K = Q(√−d) has discriminant D = −d if
d ≡ 3 (mod 4) and D = −4d otherwise. Its class number h(K) equals the number of
SL₂(ℤ)-classes of primitive positive-definite forms ax² + bxy + cy² with
b² − 4ac = D.

Every class contains exactly one **reduced** form, meaning |b| ≤ a ≤ c, with
b ≥ 0 if |b| = a or a = c. Reduced forms satisfy 3a² ≤ |D|, so counting them
is a finite search (`class_number.reduced_forms`). For example:

| d | D | h | reduced forms |
|---|---|---|---|
| 5 | −20 | 2 | (1,0,5), (2,2,3) |
| 23 | −23 | 3 | (1,1,6), (2,±1,3) |
| 47 | −47 | 5 | (1,1,12), (2,±1,6), (3,±1,4) |
| 163 | −163 | 1 | (1,1,41) |

As an independent check, `class_number_analytic` uses Dirichlet's class number
formula h(D) = −(w / 2|D|) Σ_{a=1}^{|D|} (D/a)·a. The tests confirm both
methods agree for every squarefree d < 600.

## 2. The nine Heegner numbers

Searching every squarefree d ≤ 100 000 (about 6 minutes) gives exactly

    h(Q(√−d)) = 1  ⟺  d ∈ {1, 2, 3, 7, 11, 19, 43, 67, 163}.

**A finite search can't prove the list is complete.** That's the
Heegner–Baker–Stark theorem (Heegner 1952, Baker and Stark 1967). Gauss
conjectured it; Heilbronn showed h(D) → ∞, so only finitely many exist, and
the hard part was ruling out a hypothetical tenth one. The computation only
shows there's no tenth one below the search limit.

## 3. Why e^(π√163) is almost an integer

The j-invariant is a modular function for SL₂(ℤ). With q = e^(2πiτ),

    j(τ) = E₄(τ)³ / Δ(τ) = 1/q + 744 + 196884 q + 21493760 q² + …

where E₄ = 1 + 240 Σ σ₃(n) qⁿ and Δ = q Π (1 − qⁿ)²⁴.

**CM theory.** If τ generates the ring of integers of K, then j(τ) is an
algebraic integer of degree h(K). Its conjugates are the values j at the h
reduced forms. So when h(K) = 1, j(τ) is an ordinary integer. The computation
(`j_invariant.py`, 60 digits) confirms this, and each value is a perfect cube:

| d | j |
|---|---|
| 1 | 1728 = 12³ |
| 2 | 8000 = 20³ |
| 3 | 0 |
| 7 | −3375 = −15³ |
| 11 | −32768 = −32³ |
| 19 | −884736 = −96³ |
| 43 | −884736000 = −960³ |
| 67 | −147197952000 = −5280³ |
| 163 | −262537412640768000 = −640320³ |

For contrast, d = 5 has h = 2 and j(√−5) = 632000 + 282880√5 ≈ 1264538.909,
which is not an integer.

**The near-integer.** Take τ = (1 + √−163)/2. Then q = e^(πi)·e^(−π√163) =
−e^(−π√163), which is about −3.8 × 10⁻¹⁸. Substituting:

    −640320³ = j(τ) = −e^(π√163) + 744 − 196884 e^(−π√163) + O(e^(−2π√163))

so

    e^(π√163) = 640320³ + 744 − 196884 e^(−π√163) + …
              = 262537412640768744 − 7.4992740280… × 10⁻¹³

Numerically:

    e^(π√163)  = 262537412640768743.999999999999250072597…
    640320³+744 = 262537412640768744
    difference  = −7.49927402801814 × 10⁻¹³

The gap is below 10⁻¹² because two things line up. First, h = 1 makes j(τ)
an exact integer. Second, 163 is large, so q is tiny and the 196884 q term is
the only correction that matters. The same thing happens more weakly for 43
and 67: e^(π√43) ≈ 960³ + 744 − 0.00022 and e^(π√67) ≈ 5280³ + 744 − 0.0000013.
