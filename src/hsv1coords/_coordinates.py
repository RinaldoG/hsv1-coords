"""What this library accepts as a coordinate, in one place.

Round 17: the same three-line type guard, with the same comment, had been pasted
into 17 places across four modules, so the predicate was wrong in 17 places at once.
It was `isinstance(x, int) and not isinstance(x, bool)`, which refuses `numpy.int64`,
and the audience for this library is people reading positions out of VCFs and
coordinate tables, where a numpy integer is the ordinary type: `lift(np.int64(60000),
...)` raised where the same call with a Python int answered 59,941. Before round 16
added the missing guard to the reference-source path, that call worked.

A coordinate is any integral value that is not a bool, and it is converted to `int`
so that nothing downstream has to care which kind arrived. `bool` stays refused:
`True` is 1 and a position of 1 asked by accident is worth a refusal
(standing rule 46: a naming rule, or here a type rule, lives in one function).
"""
from __future__ import annotations

import numbers


def as_coordinate(value: object, what: str = "a coordinate") -> int:
    """`value` as an `int`, or a TypeError naming what arrived."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise TypeError(f"{what} must be an int, got {value!r}")
    return int(value)
