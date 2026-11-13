"""hsv1-coords: a coordinate and annotation layer for herpes simplex virus 1."""
from .liftover import Block, Crosswalk, lift
from .repeats import CoordinateStatus, DiploidPair, Region, RepeatMap
from .resolver import (AmbiguousName, Dictionary, Entry, RepeatPairedName, normalise,
                       repeat_pair_of)

__all__ = ["AmbiguousName", "Block", "CoordinateStatus", "Crosswalk", "Dictionary",
           "DiploidPair", "Entry", "Region", "RepeatMap", "RepeatPairedName", "lift",
           "normalise", "repeat_pair_of"]
