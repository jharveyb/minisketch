# Benchmarking

## Build

```bash
cmake -B build -DCMAKE_CXX_FLAGS="-g -O2" -DMINISKETCH_BUILD_BENCHMARK=ON
cmake --build build -j4
```

## Prerequisites

These are probably available via your package manager:

- cpupower
- cpuinfo
- util-linux

## Host setup

'Robust' benchmarking can be tricky. By default, this benchmark tool will run each combination
of parameters 10 times, and return the lowest value of measured time. You can change this
with the `--iters` flag. Besides that:

- Try to run benchmarks on bare-metal, otherwise-idle machines.
- Try to have the processor clockspeed fixed during the benchmark. There are many ways to
do this; the [cpupower](https://packages.debian.org/trixie/linux-cpupower) tooling is a good option.

```bash
# Check current settings (clockspeed, governor, etc.)
sudo cpupower -c all info

# Set governor to performance
sudo cpupower -c all -g performance

# Lock clockspeeds to 1.8 GHz (separate from changing the governor)
sudo cpupower -c all -d 1800MHz
sudo cpupower -c all -u 1800MHz

# For newer x86 machines, you may also want/need to change performance hints like
# EPP, EPB, or CPPC settings.
```

For further info: <https://wiki.archlinux.org/title/CPU_frequency_scaling>

Separate from clockspeed, your processor may have heterogenous core types. In that case,
processes may migrate between core types, which would affect benchmark results.
The [cpuinfo](https://github.com/pytorch/cpuinfo) tool can show this info. For example,
on a Rockchip RK3399:

```bash
Packages:
        0: Unknown
Microarchitectures:
        2x Cortex-A72
        4x Cortex-A53
Cores:
        0: 1 processor (0), ARM Cortex-A72
        1: 1 processor (1), ARM Cortex-A72
        2: 1 processor (2), ARM Cortex-A53
        3: 1 processor (3), ARM Cortex-A53
        4: 1 processor (4), ARM Cortex-A53
        5: 1 processor (5), ARM Cortex-A53
Clusters:
        0: 2 processors (0-1),  0: 2 cores (0-1), ARM Cortex-A72
        1: 4 processors (2-5),  1: 4 cores (2-5), ARM Cortex-A53
Logical processors (System ID):
        0 (4)
        1 (5)
        2 (0)
        3 (1)
        4 (2)
        5 (3)
```

Now we can use taskset from [util-linux](https://github.com/util-linux/util-linux) to pin the benchmark to one of the A72 cores (taskset seems to use the System ID for cores):

```bash
taskset -c 5 ./build/bin/bench --sweep --sweep-errors --data=4096
```

## Output format

`bench` writes long-format TSV to stdout with a single header row:

```
metric	bits	capacity	errors	data_len	implementation	value
```

One row per (metric, configuration, implementation) tuple. Implementations (`GENERIC`, `CLMUL`, `CLMUL_TRI`) appear as a column value rather than as separate columns, so results can be piped straight into a plotting tool. Rows for implementations that are unavailable (not compiled in, not supported by the running CPU, or not defined for a given field size) are simply omitted. Informational notes (e.g. "--errors ignored in --sweep-errors mode") go to **stderr** so stdout stays pure TSV.

Three metrics are emitted:

* `recover[ms]` — best-of-`iters` wall time for a single `minisketch_decode` call.
* `create[ns]` — best-of-`iters` wall time divided by `data_len * capacity`. `minisketch_add_uint64` XORs the element into every one of the sketch's `capacity` syndrome accumulators, so per-call cost scales linearly with capacity; this normalisation exposes the per-syndrome primitive cost and should be roughly constant across capacities for a given implementation.
* `create-total[ms]` — best-of-`iters` absolute wall time for the full `data_len` add calls; use this when you care about raw throughput rather than primitive cost.

## CLI flags

Run `bench --help` for the full list. The key flags are:

* `--syndromes N`, `--errors N`, `--iters N` — per-run configuration (sketch capacity, number of distinct elements added, benchmark iterations).
* `--impl N` — restrict to a single implementation (0 = `GENERIC`, 1 = `CLMUL`, 2 = `CLMUL_TRI`). Default iterates every implementation compiled into the library.
* `--data N` — override the number of add calls in the `create` phase (default: `errors * 10`).
* `--sweep` — replace the default `bits = 2..64` loop with a fixed matrix of `bits ∈ {64}` × `capacity ∈ {1024, 2048, ..., 8192}`.
* `--sweep-errors` — implies `--sweep`; varies `errors` over values from 8 to 2048.
* `--create-only` — skip the `recover` benchmark.
* `--data-sweep` — implies `--create-only`; varies `data_len` over a fixed list of sizes instead of using `--data`.

### Examples

Baseline single-configuration run (defaults: capacity 150, errors 150, iters 10, all impls, all bit sizes from 2 to 64):

```bash
./build/bin/bench
```

Test all implementations, for bit size 64, across many combinations of capacity and error count, and save the results:

```bash
./build/bin/bench --sweep --sweep-errors --bits 64 --data=1024 > bench_64bit_results.tsv
```

Pin to the generic implementation for a larger capacity, comparing only one TSV row per bit size:

```bash
./build/bin/bench --impl=0 --syndromes=4096 --errors=1024 --iters=5
```

Fixed-matrix sweep, three iterations per point:

```bash
./build/bin/bench --sweep --iters=3
```

Also sweep error counts within each capacity (useful for studying how decode time scales with sketch occupancy):

```bash
./build/bin/bench --sweep --sweep-errors --iters=3
```

Time only the add path (skip decode) with a specific number of adds per sketch:

```bash
./build/bin/bench --sweep --create-only --data=10000 --iters=5
```

Sweep the number of adds to probe the add-path at different data volumes:

```bash
./build/bin/bench --sweep --data-sweep --iters=3
```

### Comparing CLMUL and non-CLMUL hosts

CLMUL support is compiled conditionally (see [`cmake/SystemIntrospection.cmake`](cmake/SystemIntrospection.cmake)) and is gated at runtime via a CPUID check in `src/minisketch.cpp`. To compare a CLMUL-capable host against one without CLMUL, build the binary once on the CLMUL host and run the same binary on both:

* On a CLMUL host, `GENERIC`, `CLMUL`, and (where applicable) `CLMUL_TRI` rows are emitted.
* On a non-CLMUL host, `minisketch_create` returns NULL for the CLMUL impls and only `GENERIC` rows appear.
* On an ARM host, `CLMUL` represents hosts with PMULL instruction support.

`CLMUL_TRI` is only available for field sizes that admit an irreducible trinomial over GF(2); 32-bit and 64-bit fields do not, so the `--sweep` matrix never produces `CLMUL_TRI` rows. To benchmark `CLMUL_TRI`, pick a nearby trinomial-friendly size (for example `--impl=2 --syndromes=4096 --errors=4096` with bit sizes 31 or 63 via the default non-sweep loop).

## Plotting

Read the end of the `doc/plot_bench.py` script for flag definitions. To set up dependencies, if you don't have them system-wide:

`uv sync`

Dependencies are just matplotlib + pandas. These assume you're running commands from the `doc/` directory.

Plot the decode time of multiple sketch difference counts, for certain sketch capacities, for all implementations:

```bash
uv run plot_bench.py --plots=diff --fixed-capacity 128,512,2048 $BENCH_RESULTS
```

Plot the decode time for multiple sketch capacities, where 50% of the sketch is differences, for the 58 bit field:

```bash
uv run plot_bench.py --plots=cap50 --bits=58 $BENCH_RESULTS
```

Plot the decode time for multiple (sketch difference count, sketch capacity) pairs, with capacity as the x-axis:

```bash
uv run plot_bench.py --fixed-errors 16,64,256,512,1024 --plots=cap-fixed $BENCH_RESULTS
```

Compute the speedup of the CLMUL implementation vs. the GENERIC one over all (difference count, capacity) combinations:

```bash
uv run plot_bench.py --plots=stats $BENCH_RESULTS
```
