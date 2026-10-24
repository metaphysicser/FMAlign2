"""Exercise symbol handling with synthetic data and a controlled backend.

Usage: python tests/test_symbols_cli.py /path/to/fmalign2
       python tests/test_symbols_cli.py /path/to/symbols_pipeline_test --pipeline
The backend tests record/gap/content handling, not alignment quality.
"""
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
pipeline = "--pipeline" in sys.argv[2:]
mock = Path(__file__).with_name("mock_aligner.py").resolve()
rng = random.Random(20261024)
ancestor = "".join(rng.choice("ACGT") for _ in range(500))
symbols = "RYSWKMBDHVN"
records = []
for i in range(4):
    sequence = list(ancestor + "ACGT" * i)
    for j, char in enumerate(symbols):
        sequence[j] = char
        sequence[240 + j] = char
        sequence[-1 - j] = char
    sequence[100] = "U"
    # Duplicate user IDs are safe because fragments receive unique internal IDs.
    records.append((f"sample description {i}", "".join(sequence)))

def read_fasta(path):
    result = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            result.append([line[1:], ""])
        elif line.strip():
            result[-1][1] += line.strip()
    return result

results = []
with tempfile.TemporaryDirectory(prefix="fmalign2-symbols-") as directory:
    root = Path(directory)
    for mode in ("reorder", "anchor", "drop", "duplicate", "corrupt", "truncate"):
        case = root / mode
        case.mkdir()
        (case / "input.fasta").write_text("".join(f">{name}\n{seq.lower()}\n" for name, seq in records))
        prefix = "call " if os.name == "nt" else ""
        template = prefix + f'"{sys.executable}" "{mock}" {{input}} {{output}} --mode {mode}'
        (case / "command.txt").write_text(template)
        command = [str(binary), "-i", "input.fasta", "-o", "output.fasta", "-p", "command.txt",
                   "-w", "work", "-t", "2", "-l", "2000"]
        if pipeline:
            command = [str(binary), "input.fasta", "output.fasta", "command.txt", mode]
        run = subprocess.run(command, cwd=case, capture_output=True, text=True, timeout=60)
        if mode in ("reorder", "anchor"):
            assert run.returncode == 0, (mode, run.returncode, run.stdout, run.stderr)
            output = read_fasta(case / "output.fasta")
            assert len({len(seq) for _, seq in output}) == 1
            assert [(name, seq.replace("-", "")) for name, seq in output] == records
            assert any("-" in seq for _, seq in output), "Fixture did not exercise gaps"
        else:
            assert run.returncode != 0, (mode, "Corrupt output accepted")
            assert "invalid external alignment" in run.stderr, (mode, run.returncode, run.stderr)
            assert not (case / "output.fasta").exists(), mode
        results.append({"case": mode, "exit": run.returncode, "check": "passed"})
print(json.dumps(results, indent=2))
