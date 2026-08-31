# CSO phase-state and queueing simulation revision

Updated 2026-08-31. The previous project, including its data and result files,
is saved in `HarqDelayPrediction-main_before_CSO_20260831.zip`.

## Scope

This revision implements the manuscript's auxiliary CSP-phase DTMC and
aligns the packet queueing simulation with decoding-completion release,
FIFO retransmission priority, and a complete generated-packet cohort.
The waiting and service formulas remain phase-averaged approximations;
they are not replaced by a per-process analytical model.
Existing figures and CSV results have not been recalculated or overwritten.

The source files before this numerical change are saved separately in
`HarqDelayPrediction_before_phase_fifo_20260831.zip`. Scenario values are
preserved. The relevant sweep scripts now use separate output folders so
old curves are not overwritten or presented as outputs of the new code.

## File roles and states

| File | Role |
|---|---|
| `getDvpHarqCSO.m` | Proposed model: combine auxiliary overflow, waiting, and service in the total DVP. |
| `getQLenPmfCSO.m` | Proposed model: stationary auxiliary queue PMF from the `(q,m,k_xi)` DTMC. |
| `getWaitPmfCSO.m` | Proposed model: waiting PMF given admission, with a phase-averaged blockage weight. |
| `getServiceDvpCSO.m` | Proposed model: phase-averaged service violation using feasible-attempt counts. |
| `simulate_queue_IR_cso.m` | Queueing simulation: sample packet arrivals and outcomes, and process periodic blocked slots. |
| `blerVsNrb_vsDelayTarget.m`, `Rev_blerVsNrb_vsDelayTarget.m` | Proposed-model sweeps; save `*Theory.csv`. |
| `blerVsNrb_vsDelayTarget_queuingsim.m`, `Rev_blerVsNrb_vsDelayTarget_queuingsim.m` | Queueing-simulation sweeps; save `*Simulation.csv`. |

In `getQLenPmfCSO`, the MATLAB index is
`s=k_xi*(Qmax+1)*M+q*M+m`, one greater than the paper's zero-based index.
`q` is the auxiliary maintained-packet count and `m` selects the conditional
PER. `CSO=1/xi` specifies an integer period. Phase zero is blocked; every
transition advances the phase modulo `xi`. A stored `(0,m>1,k_xi)` row goes
to `(0,1,next phase)` and has zero stationary probability. `CSO=0` uses one
data phase without blockage.

The solver forms the transition over one CSP cycle, solves its stationary
distribution, and propagates that distribution through the individual
phases. It then averages the joint probabilities over `m` and phase to
obtain `QLenPmf`. It does not replace the phase-dependent matrix by an
average of the blocked and data-slot transitions. For endpoint parameters
with several closed classes, their weights are determined from the empty
state. An optional second output provides the sparse full matrix, joint
stationary PMF, valid indices, and stationary residual for inspection.

In `simulate_queue_IR_cso`, `k_xi` is the CSP **length**, not the phase.
The current phase is `mod(slotId,k_xi)`, with zero selecting the blocked
branch. `slotId=1` is the first simulated slot, preserving the previous
initial phase. Ready packet IDs are stored in `readyIds`, attempt counts
in `TxCount`, and pending decoding/feedback events in ring buffers. The
simulation uses the manuscript's nonblocking HARQ-process assumption;
no separate finite-H process-allocation rule is added.

At each slot boundary, the simulator performs these operations in order:

1. Apply decoding outcomes whose processing interval has ended. Release
   population capacity on success or final-attempt failure.
2. Make packets whose feedback has arrived ready for retransmission.
3. Generate arrivals, while within the configured arrival horizon, and
   reject an arrival if the maintained population is full at that boundary.
4. In a data slot, schedule at most one ready packet. Pending retries have
   priority and are FIFO by feedback-availability time. First transmissions
   are FIFO by arrival. With the priority flag disabled, all ready packets
   use one FIFO list by readiness time.

Transmission in slot `s` completes decoding at the start of slot
`s+1+ceil(decodingDelay)`. A failed attempt below the cap can retransmit
from slot `s+1+ceil(decodingDelay)+ceil(feedbackDelay)`. A NACK available at
that boundary can be used immediately. Outcomes are sampled at transmission
but do not affect capacity or scheduling until their decoding event.
CSO blocks data transmission, not decoding completion or feedback.

The simulation's waiting output is arrival-to-first-transmission time, not
the auxiliary waiting PMF from `getWaitPmfCSO`. Its DVP scripts use total
latency and admission drops. Discards have infinite latency. `simLen`
specifies the number of slots in which arrivals are generated, not a fixed
number of packets. After that horizon the loop drains all admitted packets
to successful decoding or discard without adding arrivals. Therefore
`numel(totalLatency)+NumDroppedPkts` counts all generated packets. The
queue-size output retains only the original `MaxNumSlots` observations.
The optional sixth output contains packet/event records and cohort counts.
The existing nine-input call keeps one decoding slot; an optional tenth
input allows a separate raw decoding delay. Both stage delays are rounded
up separately, as in the service-delay calculation.

### Boundary and measurement conventions

At the full auxiliary buffer, an arriving packet is rejected even if a
packet leaves in the same modeled slot. A success sends `(Qmax,m)` to
`(Qmax-1,1)` without an arrival factor; at `m=M`, either outcome does so.
The packet simulator releases completed packets at the boundary before
testing a new arrival, but a transmission later in that slot cannot free
capacity for the arrival. The DTMC still omits individual decoding and
feedback timers, as required by the auxiliary approximation.

For `M=1`, the MATLAB queue function uses only `m=1`. From the empty state,
the data-slot self-transition has probability one; in a blocked slot,
the probabilities to `(0,1)` and `(1,1)` are `1-f` and `f`. Nonempty states
use the final-attempt branch. This follows the corrected Table 2 rows:
there is no attempt-2 state and the non-final first-attempt self row is
not applied again when `M=1`.

Terminal discards and overflow drops violate every finite delay target.
For delivered packets, the sweep scripts use `totalLatency > d/T`, so
completion at the deadline is timely and the final ACK is excluded.
The separate 5G postprocessor, including its additional 0.5-ms allowance,
was not changed by this MATLAB revision. Its timestamp convention and
observation horizon need to be checked independently against the metric.

No run-level confidence interval or approximation-error summary is
generated by the queueing function. Such results require the corresponding
run-level estimates and matched analysis/simulation outputs. The formula
tests listed below do not establish those empirical results.

## Waiting delay

- `getWaitPmfCSO` returns the PMF conditional on admission. Entry `k+1`
  corresponds to `k` waiting slots; the PMF sums to one.
- `pDrop` remains the auxiliary model's `pi_Qmax`, used as an overflow
  approximation rather than a physical buffer-occupancy identity.
- The one-slot shift uses `(1-beta_k)*fu(k) + beta_(k-1)*fu(k-1)`.
  An additional array entry retains the last shifted probability.
- Per-packet attempt counts are combined by convolution. This implements
  the waiting-count distribution without stopping at an intermediate zero.
  At the attempt cap, both success and discard remove a packet. This also
  handles `M=1` without changing the shared recursive waiting function.
- The output length is `M*(Qmax-1)+2`, including the zero-wait entry.
  If admission has zero stationary probability, the conditional PMF is
  undefined and the function reports an error.

## Service delay and total DVP

The existing eight-argument call remains supported:

```matlab
dvp = getDvpHarqCSO(d, delta, pHarq, f, T, M, Qmax, CSO);
```

`delta` is feedback delay in slots. The omitted decoding delay defaults to
one slot, as in the previous `delta+2` RTT expression. To pass it separately:

```matlab
zeta = 0.5;  % Raw decoding delay in slots
delta = 0.5; % Raw feedback delay in slots
dvp = getDvpHarqCSO(d, delta, pHarq, f, T, M, Qmax, CSO, zeta);
```

`getServiceDvpCSO` applies `ceil` separately to decoding and feedback delays,
then uses `RTT = 1 + tau_d + tau_f`. Thus the example gives an RTT of three
slots. Do not add `zeta` to `delta` or subtract it from `d` at the call site.
`d` and `T` must share a time unit; the stage delays are already in slots.

For residual budget `x`, the two attempt counts are:

```matlab
M_noCS = min(M, max(0, floor((x + tau_f) / RTT)));
M_CS   = min(M, max(0, floor((x + tau_f - 1) / RTT)));
alpha  = min(1, M_noCS * CSO);
```

An empty failure product is one. The service violation is the mixture of
the two failure products with weights `1-alpha` and `alpha`. The final
successful decoding need not wait for its feedback to meet the deadline.
Failure through all `M` attempts remains a violation even for a long deadline.

The total DVP includes the auxiliary drop probability, the admitted waiting
tail beyond the deadline, and the remaining service violation. The last two
terms are multiplied by the admission probability once.

`CSO` means `1/xi`, not `xi`; use zero to disable blockage. The functions
report an error outside `M*Qmax < xi` (waiting) or `M*RTT < xi` (service),
where the single-blockage approximation is not supported.

## Script changes

| Script | Output root after this revision |
|---|---|
| `blerVsNrb_vsDelayTarget.m` | `results_phase_fifo` |
| `blerVsNrb_vsDelayTarget_queuingsim.m` | `results_phase_fifo` |
| `Rev_blerVsNrb_vsDelayTarget.m` | `results_rev_phase_fifo` |
| `Rev_blerVsNrb_vsDelayTarget_queuingsim.m` | `queueing_results_rtt3_phase_fifo` |
| `testnew.m` | `results_phase_fifo` |

The CSV filenames retain the `Theory`/`Simulation` distinction. These
scripts must be rerun to obtain curves from the revised implementations.
Plotting scripts that read the old output roots need to be pointed to the
new folders. Their source data and existing plots have not been changed.

- `blerVsNrb_vsDelayTarget.m` and `Rev_blerVsNrb_vsDelayTarget.m` now pass
  `zeta` to the service calculation.
- In the revision script's attempt-cap and queue-size sweeps, `delta` was
  changed from 4 to 1 as requested. With `zeta=1`, this gives RTT=3 slots.
  The comments reflect this setting. The separate RTT sweep still uses
  `[3 4 5 6]`.
- `testnew.m` no longer combines decoding with feedback before calling the
  analysis. Its configured feedback=1 and decoding=1 now give RTT=3.
- Korean comments and Korean error messages in eight MATLAB files were
  translated into English. Their link-level calculations and tables are
  otherwise unchanged.

## Checks

Run from the project folder in MATLAB:

```matlab
results = runtests({'testCsoDelays.m', 'testCsoQueueDtmc.m', ...
    'testCsoQueueSimulation.m'});
assert(all([results.Passed]));
```

The tests cover waiting-count enumeration, agreement with the previous
recursion for `M>1`, admission normalization, the last shifted mass, `M=1`,
zero and unit PERs, separate rounding, attempt-completion thresholds, empty
products, waiting tails, discard limits, probability ranges, monotonicity,
input shapes, time units, approximation-domain checks, and the existing
stored PER inputs. Additional tests compare every DTMC row with a separate
transcription of Table 2, compare the cycle solver with a full constrained
stationary solve, and verify phase advance, unused states, full-buffer
rejection, and sparse storage at the production CSP lengths.

Simulation tests cover decoding-time capacity release, FIFO retries across
blocked slots, immediate use of available feedback, separate rounding,
draining, arrival conservation, and comparison with an independent packet
ledger that scans and sorts events instead of using ring buffers and a
maintained priority list. These are implementation checks; they do not
replace regenerated empirical comparisons or run-level confidence intervals.

The 35 tests passed in MATLAB R2025a on 2026-08-31. A separate 100,000-slot
check using the stored PER input for 10 dB/MCS 9 generated 33,465 packets;
all were admitted, 23 were discarded, and the remaining cohort finished
with a three-slot drain. A 20,000-slot high-load check generated 15,967
packets, admitted 5,646, rejected 10,321, and drained in 15 slots. These
short checks verified cohort accounting and event handling; they were not
used to replace the paper's experimental curves or estimate rare-event CIs.
The first check used `rng(30,'twister')`, Bernoulli arrivals with `f=1/3`,
`xi=80`, `Qmax=16`, `M=4`, and one slot each for decoding and feedback.
The high-load check used seed 41, `f=0.8`, `p=[0.9,0.6,0.4,0.2]`, `xi=4`,
`Qmax=4`, `M=4`, and the same stage delays. Both enabled FIFO retry priority.
