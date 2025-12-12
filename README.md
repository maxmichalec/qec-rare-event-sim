# A Rare-Event Simulation Framework for Quantim Error Correction Codes

This project aims to implement the method for simulating and evaluating QEC codes under a circuit noise model 
as described by [Rare Event Simulation of Quantum Error-Correcting Circuits](https://arxiv.org/abs/2509.13678v2)
and use [Stim](https://github.com/quantumlib/Stim) for circuit modeling and
[PyMatching](https://github.com/oscarhiggott/PyMatching) for decoding and error correction.

The `StimCircuits/` directory is taken from [StimCircuits](https://github.com/oscarhiggott/StimCircuits)
which provides a Python implementation of some standard Stim QEC circuits, with some extensions to
support building circuits with specified deterministic faults.

### Overview
Source tree:
```text
.
├── Makefile
├── README.md
├── requirements.txt
├── include
│   ├── circuit.hh
│   ├── nlohmann
│   │   └── json.hpp
│   └── utils.hh
├── src
│   ├── circuit.cc
|   ├── rare_event_simulation.cc
│   └── python
│       ├── seed_monte_carlo.py
│       ├── sim_stim_decode.py
│       ├── sim_stim_faults.py
│       └── stim_monte_carlo.py
└── StimCircuits
    ├── ...
    ├── setup.py
    └── stimcircuits
        ├── __init__.py
        ├── surface_code.py
        └── test_surface_code.py
```

`rare_event_simulation.cc` is the main file implementing the rare-event simulation framework, including
launching the inital Monte Carlo simulation (`seed_monte_carlo.py`) which uses Stim and PyMatching,
generating the splitting sequence of consecutive physical error rates $p_i$, running the Metropolis
routing to generate a RIMC for each $p_i$, and estimating the ratio between consecutive logical error rates.

`circuit.{cc,hh}` defines the QEC circuit being simulated in terms of its gates and their possible faults.
Currently, a parameterizable distance-$d$ rotated surface code is implemented.

#### Simulation of Circuit under set of Faults
During the creation of a RIMC in the Metropolis routine, candidate *Events* of (gate, fault) pairs are sampled
and need to be checked for whether they produce a logical fault or not. This requires simulation of the QEC
circuit given a predefined set of faults to occur in Stim to obtain the resulting stabalizer measurements
and pass them to PyMatching for decoding.

When an event $E$ needs to be checked, `rare_event_simulation.cc` first generates a JSON file which describes
the set of faults that occur over the $d$ rounds of stabalizer measurements in the circuit.
```json
{
    "d": d,
    "faults: {
        "0" {
            "PREP": {
                "<q>": "<f>",
            },
            "CNOT": [
                {
                    "<q>": "<f>",
                } | null,
                {
                    "<q>": "<f>",
                } | null,
                {
                    "<q>": "<f>",
                } | null,
                {
                    "<q>": "<f>",
                } | null
            ],
            "MEASURE": {
                "<q>": "<f>",
            }
        }
        ...
        "<d-1>": {...}
    }
}
```
Each fault is identified by qubit id `q` and fault type `f`, e.g. `"X"` for a single-qubit gate.

Under the implemented circuit noise model, faults can occur from:
- `"PREP"`: state preparation faults on stabalizer measurement qubits
- `"CNOT"`: 2-qubit faults during CNOT gate operation between stabalizer and data qubits
- `"MEASURE"`: basis flips during measurement of stabalizer qubits

After this file is created, the `sim_stim_faults.py` script is launched and given the path to the faults file.
This script reads the set of faults, builds the circuit with the associated error channels inserted and 
occurring with a probability of 1. `StimCircuits/stimcircuits/surface_code.py` was extended to implement
this. A circuit model of the QEC circuit being simulated with all possible error 
channels is also created, from which a `stim.DetectorErrorModel` is extracted and is used by PyMatching to
build the underlying error graph for decoding.

The circuit with the inserted faults is then simulated by Stim, generating the stabalizer measurements which
are passed to an instance of the PyMatching MWPM decoder. The script returns whether $E$ caused a logical
error or not.

### Building
To set up the environment to run simulations, first create a Python virtual environment in the project root
directory and activate it:
```bash
python3 -m venv .venv
source .venv/bin/activate
```

Next, install the required Python packages specified in the requirements file:
```bash
pip install -r requirements.txt
```
To install this custom version of `StimCircuits`:
```bash
cd StimCircuits/
pip install -e .
```

To build the rare-event simulation executable, from the project root run:
```bash
make
```
This will generate the executable `bin/rareEventSim`.


### Running Simulations
```bash
./bin/rareEventSim d p0 pTarget [N]
```
Where `d` is the code distance, `p0` is the initial physical error rate to generate an initial logical error
rate estimation from, `pTarget` is the target logical error rate to be estimated, and `N` (optional) is the 
number of failing events sampled per iteration of the Metropolis routine.

### Current Functional Correctness Bug
The current source code implements the end-to-end simulation, but there is an implementation issue likely
with respect to the Metropolis routine or ratio estimation, so logical error estimates are produced but
are not currently accurate.
