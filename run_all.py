"""
Reproduce every result and figure from a clean checkout.

    pip install -r requirements.txt
    python run_all.py

Runs, in dependency order:
  1. energy_budget.py        -> energy_budget.csv
  2. radiotrophic_model.py   -> FBA experiment CSVs + radiotrophic_cell_model.json
  3. kinetic_model.py        -> k1..k4 CSVs
  4. dsup_analysis.py        -> dsup_* CSVs + JSON
  5. make_figures.py         -> figures/*.png
"""

import subprocess
import sys

STEPS = [
    "energy_budget.py",
    "radiotrophic_model.py",
    "kinetic_model.py",
    "dsup_analysis.py",
    "make_figures.py",
]


def main():
    for script in STEPS:
        print(f"\n{'='*64}\n>>> {script}\n{'='*64}", flush=True)
        result = subprocess.run([sys.executable, script])
        if result.returncode != 0:
            print(f"\nFAILED: {script} (exit {result.returncode})")
            sys.exit(result.returncode)
    print("\nAll results and figures regenerated.")


if __name__ == "__main__":
    main()
