
## Reference publication
https://arxiv.org/abs/2502.08789 

## Overview
MATLAB scripts to study Delay Violation Probability (DVP) for ARQ/HARQ under finite-blocklength assumptions. 
The top-level script is A0_Script_runAllCodes.m. It resets figures and runs a sequence of all plot generators based on what is commented and what is not. (see **Reproducing figures**).
This script could be used as a starting point to refer to all other script and the internal functions contained within. The name of the subscripts is mostly self-explanatory after going through the README.

### CSO model and simulation

The proposed model is implemented by `getDvpHarqCSO`, `getQLenPmfCSO`,
`getWaitFunctionCSO`, and `getServiceDvpCSO`. The model calculates the DVP
upper bound from the packet-count distribution, waiting-delay function, and
service-delay result. It uses the state `(q,m,k_xi)`, where `q` is the
aggregate packet count, `m` selects the attempt-dependent PER, and `k_xi`
tracks the CSP phase. `Qmax` is the maximum aggregate packet count.

The MATLAB Queueing Simulation is implemented separately by
`simulate_queue_IR_cso`. It allows another packet to be transmitted while a
packet is awaiting decoding or feedback, gives ready retransmissions priority,
and schedules ready retransmissions in FIFO order. Its waiting-delay output
includes contention before the initial transmission and between subsequent
transmission attempts. Its service-delay output includes the tagged packet's
transmissions, decoding intervals, and feedback intervals except after the
final attempt.

The CSO figure and sensitivity-study scripts define their respective scenario
settings. The parameter values below are general examples for the wider code
collection.

## Requirements
- MATLAB (tested on recent 2024b)
- Might not require any specialised toolboxes except porobably the basic ones like communication toolbox for some in-built functions. This would throw error whenever encountered.
- All `.m` files in the same folder (or on MATLAB path)

## Quick start example
```matlab
close all;
h1 = figure(1); h2 = figure(2); h3 = figure(3);
clf(h1,'reset'); clf(h2,'reset'); clf(h3,'reset');
bler_vs_dByRttRatio;

% ...uncomment the blocks you need...
```

## Parameters

| Name | Meaning | Units | Typical Default / Example |
|---|---|---:|---|
| `d` | End-to-end delay target. Suggested to NOT keep it a whole number of slot to avoid >=< confusions. | ms | e.g., `d = 5` |
| `delta` | Feedback (ACK/NACK) delay | slots | e.g., `delta = 4` |
| `zeta` | Decoding/processing delay | slots | e.g., `zeta = 1` |
| `Qmax` | Max queue size for **ARQ** | pkts | e.g., `Qmax = 100` |
| `M` | Max retransmissions for **ARQ** (incl. first tx) | count | e.g., `M = 3` |
| `Qmax_Harq` | Maximum aggregate packet count for **HARQ** | pkts | e.g., `Qmax_Harq = 16` |
| `M_Harq` | Max retransmissions for **HARQ** (incl. first tx) | count | e.g., `M_Harq = 3` |
| `T` | Slot length (as per numerology) | s | `T = 1e-3` |
| `hvar` | Channel power gain variance | – | `hvar = 1` |
| `V` | Dispersion coefficient (finite-blocklength) | – | `V = 1` |
| `SNRdb_avgRcvd` | Avg received SNR | dB | e.g., `SNRdb_avgRcvd = 5` |
| `f` | Bernoulli arrival prob. per slot (= arrival rate per slot) | – | e.g., `f = 0.2` |
| `pktLen` | Payload length. Choose so as to avoid segmentation. Refer paper. | bits | `pktLen = 8*100` |
| `Nrb` | Number of RBs used for data | RB | `Nrb = 10` |
| `Nre` | REs per slot used for **data** (= `Nrb * 180`, see note) | RE | `Nre = Nrb*180` |
| `B` | Effective occupied bandwidth (= `180e3 * Nrb`) | Hz | `B = 180e3*Nrb` |

### Notes on key parameters/assumptions
- **Sloting:** `T = 1e-3` assumes µ=0 (15 kHz SCS). Adjust if you change numerology. Other numerologies would be handled by adjusting T; however is not tested.
- **`Nre` (REs):** Uses **15 symbols** in time (excludes CP and typical control/DMRS overhead). This can exceed the usable RE count; adjust to your DMRS/PDCCH/PTRS settings if needed.  
- **`B` (Bandwidth):** Assumes full 180 kHz per RB per slot (12*15=180) is occupied by the data transmission. This is **effective used BW**, not carrier bandwidth.  
- **MCS:** Table 1 by default; see `parseMcs.m` for mapping.

## Reproducing figures
Each function generates a set of plots. The main script cleans figures between runs.

1. **DVP vs d/RTT ratio for different RTT**
   ```matlab
   clf(h1,'reset'); clf(h2,'reset'); clf(h3,'reset');
   bler_vs_dByRttRatio
   ```

2. **DVP vs Nrb for different delay targets**
   ```matlab
   clf(h1,'reset'); clf(h2,'reset'); clf(h3,'reset');
   blerVsNrb_vsDelayTarget
   ```

3. **DVP vs Nrb for different allocated resources per packet length**
   

4. **DVP vs Nrb for different pktLen**

5. **DVP vs Nrb for different RTT**

6. **DVP vs Nrb for different SNR**

7. **DVP vs PktLen for different SNR**

8. **DVP vs RTT for different SNR**
   
9. **Throughput vs f for different Nrb**

10. **Throughput vs f for different RTT**

10. **performance comparison**
For comparing performance with the baselines IF scheme and also comparing between event-based simulator and closed forms. Other figures mainly plot the results.
   ```matlab
   clf(h1,'reset'); clf(h2,'reset'); clf(h3,'reset');
   perf_comparison
   ```

## Tips
- Keep `delta` and `zeta` aligned with your PHY/stack model; they directly eat into the latency budget `d`. Even if decoding delay is say half a slot but the feedback for the failed message or the forwarding of the successful message happens at slot boundary, then provide decoding delay as full slot. Adjust accordingly for other scenarios. Similarly for feedback delay, include if there is a wait before retransmission. As part of feedback as they happens always as a pair, it is not necessary to separate them. that being said one can do so and adjust easily if they wish to.
- If modeling HARQ combining, ensure `M_Harq` logic and effective SNR accumulation are consistent across BLER/throughput functions.
- When tightening realism, reduce `Nre` to account for DMRS/PDCCH overhead and adjust `B` if partial RBs or mini-slots are used.
