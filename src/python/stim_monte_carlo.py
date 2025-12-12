import stim
import pymatching
import numpy as np
import pandas as pd
import argparse
import sinter
import matplotlib.pyplot as plt

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
    parser.add_argument("shots", type=int)
    parser.add_argument("-r", "--rounds", type=int, default=None)
    args = parser.parse_args()

    if args.rounds is None:
        rounds = args.d
    else:
        rounds = args.rounds

    surface_code_circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        rounds=rounds,
        distance=args.d,
        after_clifford_depolarization=args.p,
        # before_round_data_depolarization=args.p,
        # before_measure_flip_probability=args.p,
        # after_reset_flip_probability=args.p
    )
    # print(str(surface_code_circuit))

    n_errors = count_logical_errors(surface_code_circuit, args.shots)

    p_error_logical = n_errors / args.shots

    # print(f"Z-basis logical error for d={args.d}, p_phy={args.p}, shots={args.shots}: {p_error_logical}")
    print(f"{p_error_logical:.5e}")

def mc_sampling():
    tasks = [
        sinter.Task(
            circuit=stim.Circuit.generated(
                "surface_code:rotated_memory_z",
                rounds=d,
                distance=d,
                after_clifford_depolarization=noise,
                before_round_data_depolarization=noise,
                after_reset_flip_probability=noise
            ),
            json_metadata={'d': d, 'p': noise},
        )
        for d in [3, 5, 7, 9]
        # for noise in [0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1]
        for noise in [0.00001, 0.00005, 0.0001, 0.0005, 0.001]
    ]

    collected_stats: list[sinter.TaskStats] = sinter.collect(
        num_workers=16,
        tasks=tasks,
        decoders=['pymatching'],
        max_shots=100_000_000,
        max_errors=1000,
    )

    data = []
    for ts in collected_stats:
        data.append({
            "d": ts.json_metadata['d'],
            "p": ts.json_metadata['p'],
            "pL": ts.errors / ts.shots,
            "time": ts.seconds
        })
    df = pd.DataFrame(data)
    print(df)

    fig, ax = plt.subplots(1, 1)
    sinter.plot_error_rate(
        ax=ax,
        stats=collected_stats,
        x_func=lambda stats: stats.json_metadata['p'],
        group_func=lambda stats: stats.json_metadata['d'],
    )
    # ax.set_ylim(5e-7, 1e-0)
    # ax.set_xlim(5e-4, 1e-1)
    ax.loglog()
    ax.set_title("Rotated Surface Code Error Rates (Circuit Noise Model)")
    ax.set_xlabel("Phyical Error Rate")
    ax.set_ylabel("Logical Error Rate")
    ax.grid(which='major')
    ax.grid(which='minor')
    ax.legend()
    fig.set_dpi(120)  # Show it bigger
    plt.savefig("surface_code_error_rate.svg", format="svg")


if __name__ == '__main__':
    # main()
    mc_sampling()
