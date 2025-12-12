import stim
import stimcircuits
import pymatching
import numpy as np
import json
import argparse

'''
Gate faults list format
List of length rounds with elements:
{
  "PREP": { int measure_qubit_id: char pauli }
  "CNOT": [{ int measure_qubit_id: char pauli[2] }]  # length=4
  "MEASURE": { int measure_qubit_id: char pauli }
}
'''

def parse_fault_file(fpath):
    with open(fpath, 'r') as f:
        data = json.load(f)
    
    d = int(data["d"])

    assert("faults" in data)

    faults = []
    for r in range(d):
        faults.append({})

    for r in range(d):
        if "PREP" not in data["faults"][str(r)]:
            faults[r]["PREP"] = {}
        else:
            faults[r]["PREP"] = data["faults"][str(r)]["PREP"]

        if "MEASURE" not in data["faults"][str(r)]:
            faults[r]["MEASURE"] = {}
        else:
            faults[r]["MEASURE"] = data["faults"][str(r)]["MEASURE"]

        if "CNOT" not in data["faults"][str(r)]:
            faults[r]["CNOT"] = [None] * 4
        else:
            faults[r]["CNOT"] = data["faults"][str(r)]["CNOT"]

    # print(faults)

    return d, faults

def causes_logical_error(circt: stim.Circuit, circt_with_failures: stim.Circuit):
    sampler = circt_with_failures.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(1, separate_observables=True)

    # print(detection_events)
    # print(observable_flips)

    detector_error_model = circt.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    predictions = matcher.decode_batch(detection_events)

    if not np.array_equal(observable_flips[0], predictions[0]):
        return True
    else:
        return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fault_file", type=str)
    args = parser.parse_args()

    d, faults = parse_fault_file(args.fault_file)

    circt = stimcircuits.generate_circuit(
        "surface_code:rotated_memory_z",
        rounds=d,
        distance=d,
        after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001,
        after_reset_flip_probability=0.001
    )
    # print(circt)

    circt_with_failures = stimcircuits.generate_circuit(
        "surface_code:rotated_memory_z",
        rounds=d,
        distance=d,
        gate_faults=faults
    )
    # print(circt_with_failures)

    # diagram = circt_with_failures.diagram("timeline-svg")
    # with open("diagram.svg", 'w') as f:
    #     f.write(str(diagram))

    logical_error = causes_logical_error(circt, circt_with_failures)
    print(int(logical_error))


if __name__ == '__main__':
    main()
