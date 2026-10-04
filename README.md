# ZetaMine benchmark

This repository contains the benchmark runner for the ZetaMine frequent itemset mining implementation. The script compiles the C++ engine, downloads benchmark datasets if needed, runs the comparison baselines, and writes timing/output artifacts into the local `zm_bench` folder.

## Requirements

Before running the benchmark, make sure you have:

- Python 3.10+
- g++ on your PATH
- Optional: Java if you want the SPMF baselines enabled
- Python packages: `numpy`, `pandas`, and `matplotlib`

On Linux/macOS, install the essentials with:

```bash
python3 -m pip install numpy pandas matplotlib
```

On Windows, install a C++ toolchain such as MSYS2 and ensure `g++` is available on PATH.

## Run the benchmark

From the project root:

```bash
cd /workspaces/zetamine
python3 benchmark.py
```

The script will:

1. compile the C++ benchmark binary into `zm_bench/zm`
2. install `pyfim` if needed
3. download benchmark datasets into `zm_bench/data`
4. execute the benchmark and save results under `zm_bench/`

## Recommended quick run

To keep the run shorter, use the built-in smoke preset:

```bash
cd /workspaces/zetamine
ZM_PRESET=smoke python3 benchmark.py
```

Other presets available in the script:

- `smoke`
- `mid`
- `quick`
- `full` (default)

Example:

```bash
ZM_PRESET=quick python3 benchmark.py
```

## Optional environment variables

You can tune the run without changing the code:

```bash
ZM_PRESET=smoke
ZM_REPS=3
ZM_TIMEOUT=600
ZM_ABLATIONS=1
python3 benchmark.py
```

Meaning:

- `ZM_PRESET`: dataset preset (`smoke`, `mid`, `quick`, `full`)
- `ZM_REPS`: number of benchmark repetitions
- `ZM_TIMEOUT`: per-run timeout in seconds
- `ZM_ABLATIONS`: enable or disable ablation runs (`1` or `0`)

## Outputs

The script writes benchmark results to the workspace, including:

- `zm_bench/results.csv`
- `zm_bench/environment.json`
- `zm_bench/runtime_output.png`
- `zm_bench/runtime_count.png`
- `zetamine_results.zip` in the parent directory after completion

## Notes

- The full benchmark can take a while because it compiles C++, installs dependencies, and downloads datasets.
- If `java` is available, the script can enable SPMF baselines automatically.
- If `g++` is missing, the run stops immediately with a setup error.
