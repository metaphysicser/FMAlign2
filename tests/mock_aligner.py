"""Deterministic I/O fixture, not a biological alignment algorithm."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--mode", default="reorder")
args = parser.parse_args()
records = []
for line in Path(args.input).read_text().splitlines():
    if line.startswith(">"):
        records.append([line[1:], ""])
    elif line.strip():
        records[-1][1] += line.strip()
width = max(len(seq) for _, seq in records)
records = [[name, seq.lower().replace("u", "t").ljust(width, "-")] for name, seq in records]
# The CLI runs a separate backend availability probe before the real input.
if records and records[0][0].startswith("SEQENCE"):
    if args.mode == "drop":
        records.pop()
    elif args.mode == "duplicate":
        records[-1][0] = records[0][0]
    elif args.mode == "corrupt":
        seq = records[0][1]
        records[0][1] = ("c" if seq[0] == "a" else "a") + seq[1:]
    elif args.mode == "truncate":
        records[0][1] = records[0][1][:-1]
records.reverse()
Path(args.output).write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
