"""Exit codes and the exceptions that carry them."""
from __future__ import annotations

from enum import IntEnum


class ExitCode(IntEnum):
    OK = 0
    INTERNAL = 1            # unexpected internal error
    USAGE = 2               # usage error
    GATE = 3                # QC gate tripped (only with --fail-on)
    INPUT = 4               # input not found, unreadable or not recognised
    PARTIAL = 5             # finished, but some runs or inputs failed or were skipped
    LIBRARY = 6             # library problem: not a library, newer layout, not writable
    VALIDATION = 7          # validation failed
    RESOURCE = 8            # a resource limit was exceeded
    INTERRUPTED = 130       # interrupted; partials written so far are kept


class FlowCellarError(Exception):
    """Base class: the command line prints the message and exits with ``exit_code``."""

    exit_code: int = ExitCode.INTERNAL

    def __init__(self, message: str, *, hint: str | None = None):
        super().__init__(message)
        self.message = message
        self.hint = hint


class UsageError(FlowCellarError):
    exit_code = ExitCode.USAGE


class InputError(FlowCellarError):
    exit_code = ExitCode.INPUT


class LibraryError(FlowCellarError):
    exit_code = ExitCode.LIBRARY


class ValidationError(FlowCellarError):
    exit_code = ExitCode.VALIDATION


class ResourceLimitError(FlowCellarError):
    exit_code = ExitCode.RESOURCE


class MissingExtraError(FlowCellarError):
    """An optional dependency is needed: names the extra to install."""

    exit_code = ExitCode.INPUT

    def __init__(self, package: str, extra: str, what: str):
        super().__init__(f"{what} needs the package '{package}'.", hint=f'pip install "flowcellar[{extra}]"')


class NotBuiltError(FlowCellarError):
    """A part that is registered but is not in this build."""

    exit_code = ExitCode.INTERNAL

    def __init__(self, what: str):
        super().__init__(f"{what}: not built yet")
        self.what = what


def not_built(what: str):
    """``raise not_built("flowcellar.engine.api.quick_bundle")``"""
    return NotBuiltError(what)
