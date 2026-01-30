from pathlib import Path

seq = Path("Output.fa").read_text()
assert "GAATTC" not in seq
print("Test passed: EcoRI removed")
