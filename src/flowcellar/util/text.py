"""Text from input files, and numbers as the interface words them."""
from __future__ import annotations

import unicodedata


def clean_text(s, limit: int = 4096) -> str:
    """Untrusted text from an input file: control characters removed, NFC, truncated."""
    s = unicodedata.normalize("NFC", str(s))
    s = "".join(ch for ch in s if ch == "\n" or unicodedata.category(ch)[0] != "C")
    return s[:limit]


def fmt_count(x: float | None) -> str:
    """76.4 M; 8,270."""
    if x is None:
        return "not known"
    if abs(x) >= 1e6:
        v = x / 1e6
        return (f"{v:.0f}" if abs(v) >= 100 else f"{v:.1f}") + " M"
    return f"{int(round(x)):,}"


def fmt_bases(x: float | None) -> str:
    """186 Gb; 45.7 Gb; 540 Mb."""
    if x is None:
        return "not known"
    v = x / 1e9
    if abs(v) >= 100:
        return f"{v:.0f} Gb"
    if abs(v) >= 1:
        return f"{v:.1f} Gb"
    return f"{x / 1e6:.0f} Mb"


def fmt_pct(x: float, digits: int = 0) -> str:
    return f"{x * 100:.{digits}f} %"


def fmt_hours(h: float) -> str:
    if h < 1:
        return f"{round(h * 60)} minutes"
    return f"{h:.1f} hours" if h < 10 else f"{h:.0f} hours"
