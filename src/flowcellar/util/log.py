"""Logging: text on stderr, one line per event; the same events as JSON lines with --log-format json."""
from __future__ import annotations

import json
import logging
import sys
import time

LOGGER_NAME = "flowcellar"
log = logging.getLogger(LOGGER_NAME)


class _JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        d = {"t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(record.created)), "level": record.levelname.lower(),
             "event": getattr(record, "event", None) or "message", "message": record.getMessage()}
        d.update(getattr(record, "fields", None) or {})
        return json.dumps(d, ensure_ascii=False, default=str)


class _TextFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        msg = record.getMessage()
        if record.levelno >= logging.WARNING:
            return f"{record.levelname.lower()}: {msg}"
        return msg


def setup(verbosity: int = 0, fmt: str = "text", stream=None) -> logging.Logger:
    """verbosity: -1 quiet (warnings only), 0 normal, 1 and more debug. Safe to call again."""
    level = logging.WARNING if verbosity < 0 else logging.INFO if verbosity == 0 else logging.DEBUG
    handler = logging.StreamHandler(stream or sys.stderr)
    handler.setFormatter(_JsonFormatter() if fmt == "json" else _TextFormatter())
    log.handlers[:] = [handler]
    log.setLevel(level)
    log.propagate = False
    return log


def event(name: str, message: str, *, level: int = logging.INFO, **fields) -> None:
    """One decision or one result: ``event("run.quick", "quick report written", run_key=..., seconds=8.1)``."""
    log.log(level, message, extra={"event": name, "fields": fields})


def is_tty() -> bool:
    """No progress display when stderr is not a terminal."""
    try:
        return sys.stderr.isatty()
    except (AttributeError, ValueError):
        return False
