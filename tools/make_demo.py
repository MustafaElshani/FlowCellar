#!/usr/bin/env python3
"""Write the FlowCellar example (an invented laboratory: run folders, read data, analysis folders). Invented names only.

    python tools/make_demo.py OUT_DIR [--scale 0.01] [--jobs 8]

The same as ``python -m flowcellar.testing.demo``; see that module for what is written.
"""
import sys

from flowcellar.testing.demo import main

if __name__ == "__main__":
    sys.exit(main())
