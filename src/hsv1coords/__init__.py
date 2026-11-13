"""hsv1-coords: a coordinate and annotation layer for herpes simplex virus 1."""
from .liftover import Block, Crosswalk, lift
from .resolver import (AmbiguousName, Dictionary, Entry, RepeatPairedName, normalise,
                       repeat_pair_of)

__all__ = ["AmbiguousName", "Block", "Crosswalk", "Dictionary", "Entry",
           "RepeatPairedName", "lift", "normalise", "repeat_pair_of"]
