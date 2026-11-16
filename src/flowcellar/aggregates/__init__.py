"""The aggregate kinds of a bundle: building them, merging them, deriving numbers from them.

Every stored number is one of nine kinds; each kind has one merge rule. Only integers are accumulated, so
merging in any order is exact.

    from flowcellar.aggregates import build, merge, derive
"""
from flowcellar.aggregates import build, derive, merge  # noqa: F401
from flowcellar.aggregates.merge import NotMergeable, merge_cell, merge_cells, merge_kind, rollup  # noqa: F401

KINDS = ("scalars", "hists", "hists2d", "series", "counts", "channels", "tables", "curves", "grids")
