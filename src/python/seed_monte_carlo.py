import stim
import pymatching
import numpy as np
import argparse

# From: https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb
def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)

    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("d", type=int)
    parser.add_argument("p", type=float)
    parser.add_argument("-e", "--errors", type=int, default=1000)
    parser.add_argument("-s", "--shots", type=int, default=1_000_000)
    args = parser.parse_args()

    surface_code_circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=args.d,
        distance=args.d,
        after_clifford_depolarization=args.p,
        # before_round_data_depolarization=args.p,
        before_measure_flip_probability=args.p,
        after_reset_flip_probability=args.p
    )

    n_errors = count_logical_errors(surface_code_circuit, args.shots)

    p_error_logical = n_errors / args.shots

    print(f"{p_error_logical:.5e}")


if __name__ == '__main__':
    main()
