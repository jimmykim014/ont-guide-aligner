"""Tests for functions in the Rust core that the CLI calls directly."""

import itertools
import random

from ont_guide_aligner import _core


def brute_force_close_count(seqs):
    return sum(
        1
        for a in seqs
        if any(b != a and sum(x != y for x, y in zip(a, b)) == 1 for b in seqs)
    )


def test_close_guide_count_none():
    assert _core.close_guide_count(["AAAA", "CCCC", "GGGG", "TTTT"]) == 0


def test_close_guide_count_empty():
    assert _core.close_guide_count([]) == 0


def test_close_guide_count_one_pair():
    assert _core.close_guide_count(["ACGTACGT", "ACGAACGT", "TTTTTTTT"]) == 2


def test_close_guide_count_chain_counts_each_guide_once():
    # A-B and B-C are 1 apart, A-C are 2 apart: all three guides have a close neighbour
    assert _core.close_guide_count(["AAAA", "AACA", "GACA"]) == 3


def test_close_guide_count_matches_brute_force():
    rng = random.Random(1)
    seqs = set()
    while len(seqs) < 200:
        seqs.add("".join(rng.choice("ACGT") for _ in range(20)))
    seqs = sorted(seqs)
    # add some 1-substitution neighbours so the answer isn't trivially 0
    extra = []
    for s in seqs[:30]:
        p = rng.randrange(20)
        b = rng.choice([c for c in "ACGT" if c != s[p]])
        n = s[:p] + b + s[p + 1:]
        if n not in seqs and n not in extra:
            extra.append(n)
    seqs = seqs + extra
    rng.shuffle(seqs)
    expected = brute_force_close_count(seqs)
    assert expected >= 2 * len(extra) > 0
    assert _core.close_guide_count(seqs) == expected
