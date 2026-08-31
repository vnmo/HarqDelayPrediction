# CSO waiting/service delay revision

Updated 2026-08-31. The previous project, including its data and result files,
is saved in `HarqDelayPrediction-main_before_CSO_20260831.zip`.

## Scope

This revision updates the waiting-delay and service-delay approximations and
their combination in the total DVP. It does not change the numerical logic
of `getQLenPmfCSO.m`, the shared non-CSO waiting functions, or the queueing
simulation. Comments now identify the Proposed-model and simulation files
and describe their stored states and branches.
The queue DTMC still averages CSO in its transition probabilities; a CS-phase
state and individual HARQ-process states have not been added.
Existing figures and CSV results have not been recalculated or overwritten.

The source files before these comment changes are saved separately in
`HarqDelayPrediction_CSO_before_comments_20260831.zip`. Scenario values and
output paths were preserved during the comment changes.

## File roles and states

| File | Role |
|---|---|
| `getDvpHarqCSO.m` | Proposed model: combine auxiliary overflow, waiting, and service in the total DVP. |
| `getQLenPmfCSO.m` | Proposed model: stationary auxiliary queue PMF from the `(q,m)` DTMC. |
| `getWaitPmfCSO.m` | Proposed model: waiting PMF given admission, with a phase-averaged blockage weight. |
| `getServiceDvpCSO.m` | Proposed model: phase-averaged service violation using feasible-attempt counts. |
| `simulate_queue_IR_cso.m` | Queueing simulation: sample packet arrivals and outcomes, and process periodic blocked slots. |
| `blerVsNrb_vsDelayTarget.m`, `Rev_blerVsNrb_vsDelayTarget.m` | Proposed-model sweeps; save `*Theory.csv`. |
| `blerVsNrb_vsDelayTarget_queuingsim.m`, `Rev_blerVsNrb_vsDelayTarget_queuingsim.m` | Queueing-simulation sweeps; save `*Simulation.csv`. |

In `getQLenPmfCSO`, `s=q*M+m` indexes the auxiliary maintained-packet count
`q` and the attempt index `m` used to select the conditional PER. `CSO=1/xi`
enters the transition probabilities as a blocked-slot probability. CSP phase
is not another dimension of this matrix. The waiting and service functions
also use phase-averaged weights; they do not advance a slot clock.

In `simulate_queue_IR_cso`, `k_xi` is the CSP **length**, not the phase.
The current phase is `mod(slotId,k_xi)`, with zero selecting the blocked
branch. Ready packet IDs are stored in `txQueue_pktId`, attempt counts in
`TxCount`, and future retry availability in `reTxPktIdList`. A retry becomes
ready after transmission, decoding, and feedback delay. Other ready packets
may transmit in the meantime. No finite HARQ-process pool is allocated.

The simulation processes arrivals and retry availability before checking
whether the slot is blocked. A blocked slot prevents data transmission;
it does not undo those earlier updates. An idle branch means no packet is
ready, while a packet may still be awaiting a retry. Success or final failure
releases the admitted-packet counter when the outcome is sampled. Successful
latency includes the later decoding completion; the final ACK is excluded.

The simulation's waiting output is arrival-to-first-transmission time, not
the auxiliary waiting PMF from `getWaitPmfCSO`. Its DVP scripts use total
latency and admission drops. Discards have infinite latency; unresolved
packets are omitted from the returned latency arrays. `simLen` specifies
simulation slots, not a fixed number of generated packets.

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
results = runtests('testCsoDelays.m');
assert(all([results.Passed]));
```

The tests cover waiting-count enumeration, agreement with the previous
recursion for `M>1`, admission normalization, the last shifted mass, `M=1`,
zero and unit PERs, separate rounding, attempt-completion thresholds, empty
products, waiting tails, discard limits, probability ranges, monotonicity,
input shapes, time units, approximation-domain checks, and the existing
stored PER inputs. These are implementation checks; they do not replace
empirical validation against queueing simulation or measurements.
