#!/usr/bin/env python3
"""Write a standalone HTML file with the pocket demo (no React needed)."""

import argparse
import os
import sys

# Allow running from repo root or backend/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app.services.pocket_viz_3d import write_standalone_html  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output", default="pocket_3d_demo.html", help="Output HTML path")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    write_standalone_html(args.output, seed=args.seed)
    print(f"Wrote {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
