import stim
import pymatching
import numpy as np
import json
from enum import Enum

class Gate(Enum):
    PREP0 = 1
    PREPPLUS = 2
    CNOT = 3
    MEAS_X = 4
    MEAS_Z = 5
    IDLE = 6

class PauliOp(Enum):
    I = 1, 'I'
    X = 2, 'X'
    Y = 3, 'Y'
    Z = 4, 'Z'

def build_schedule(d):
    n_dataq = d * d
    schedule = []
    stabalizers = []

    for i in range(d - 1):
        for j in range(d - 1):
            nw = i * d + j
            ne = i * d + (j + 1)
            sw = (i + 1) * d + j
            se = (i + 1) * d + (j + 1)
            stabalizers.append([nw, ne, sw, se])

    qubit_index = n_dataq
    for r in range(d * 2):
        isXstab = (r % 2 == 0)
        for data_qubits in stabalizers:
            anc_id = qubit_index
            qubit_index += 1
            # Prepare ancilla qubit
            if isXstab:
                schedule.append({"gate": Gate.PREPPLUS, "qubits": (anc_id,)})
            else:
                schedule.append({"gate": Gate.PREP0, "qubits": (anc_id,)})

            # Stabalizer measurement with CNOTS
            for q in data_qubits:
                if isXstab:
                    schedule.append({"gate": Gate.CNOT, "qubits": (anc_id, q)})
                else:
                    schedule.append({"gate": Gate.CNOT, "qubits": (q, anc_id)})
            
            # Measurement of ancilla
            if isXstab:
                schedule.append({"gate": Gate.MEAS_X, "qubits": (anc_id,)})
            else:
                schedule.append({"gate": Gate.MEAS_Z, "qubits": (anc_id,)})
        
        # Insert idles for data qubits
        for q in range(n_dataq):
            schedule.append({"gate": Gate.IDLE, "qubits": (q,)})

    return schedule, n_dataq

def construct_and_sim(schedule, n_data, faults):
    circt = stim.Circuit()

    for op in schedule:
        gate = op["gate"]
        # fault = op["fault"]
        qubits = op["qubits"]

        # Append gate to circuit
        match gate:
            case Gate.PREP0:
                circt.append(f"RZ", qubits[0])
            case Gate.PREPPLUS:
                circt.append(f"RX", qubits[0])
            case Gate.CNOT:
                circt.append(f"CNOT", [qubits[0], qubits[1]])
            case Gate.MEAS_X:
                circt.append(f"MX", qubits[0])
            case Gate.MEAS_Z:
                circt.append(f"MZ", qubits[0])
            case Gate.IDLE:
                circt.append(f"I", qubits[0])
            case _:
                circt.append(f"// unknown gate/operation: {gate}")

        # Inject fault as Pauli gate
        # if gate is Gate.CNOT:
        #     if fault[0] is not PauliOp.I:
        #         circt.append(f"{fault[0].name}", qubits[0])
        #     if fault[1] is not PauliOp.I:
        #         circt.append(f"{fault[1].name}", qubits[1])
        # else:
        #     if fault[0] is not PauliOp.I:
        #         circt.append(f"{fault[0].name}", qubits[0])

    # sim = stim.TableauSimulator()

    # sampler = circt.compile_sampler()
    # measurements = sampler.sample(shots=1)

    return circt

def decode_pymatching(detector_error_model, measurements):
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)
    correction = matcher.decode(measurements)
    return correction

def main():
    schedule, n_dataq = build_schedule(3)
    circt = construct_and_sim(schedule, n_dataq, [])

    sampler = circt.compile_sampler()
    measurements = sampler.sample(shots=1)
    dem = circt.detector_error_model()
    print(str(circt))
    matching = pymatching.Matching.from_detector_error_model(dem)
    print(matching)
    correction = matching.decode(measurements)
    print(correction)


if __name__ == '__main__':
    main()
