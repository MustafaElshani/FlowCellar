#!/usr/bin/env python3
"""Write fake MinKNOW run folders with invented names: a thin wrapper around flowcellar.testing.fixtures.

    python tools/make_fixtures.py OUT_DIR [--reads 4000] [--hours 8] [--family minion|promethion] [--tree]
"""
import sys

from flowcellar.testing.fixtures import main

if __name__ == "__main__":
    sys.exit(main())
