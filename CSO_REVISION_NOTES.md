# CSO waiting/service delay revision

Updated 2026-08-31. The previous project, including its data and result files,
is saved in `HarqDelayPrediction-main_before_CSO_20260831.zip`.

## Scope

This revision updates the waiting-delay and service-delay approximations and
their combination in the total DVP. It does not change `getQLenPmfCSO.m`,
the shared non-CSO waiting functions, or the queueing simulation files.
The queue DTMC still averages CSO in its transition probabilities; a CS-phase
state and individual HARQ-process states have not been added.
Existing figures and CSV results have not been recalculated or overwritten.

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
