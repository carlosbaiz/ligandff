"""Per-atom charge penalty (P-II eq 4; SPEC 5.11), verified exactly against program 5.0 output (B0).

    sigma_i = sqrt( sum_n cbrt(|beta_n| + delta) * sigma_n^2 ),   delta = 0.05**6

beta_n is the net charge that one parameter instance moves onto atom i and sigma_n is that parameter's analogy
penalty; instances of parameters present in CGenFF (penalty 0) contribute nothing.
"""
import math

DELTA = 0.05 ** 6


def atom_penalty(contributions):
    """contributions: iterable of (net charge transfer onto the atom, parameter penalty)."""
    return math.sqrt(sum((abs(dq) + DELTA) ** (1 / 3) * p * p for dq, p in contributions))
