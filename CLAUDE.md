# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-page, zero-build web app that visualizes a **network coherogram** of ground/space
magnetometer data to detect ULF (Pc3–Pc5, ~2–50 mHz) waves. The app's UI, CSS, and orchestration
live in **`index.html`**; the pure DSP primitives are factored into **`dsp.js`** (a plain global
`<script>` in the browser, `require()`-able in Node) so they can be unit-tested. There is no
package.json, build step, or backend.

## Running / developing

- Open `index.html` directly in a browser, or serve it (e.g. `python3 -m http.server`) and load
  it. Serving over HTTP avoids any `file://` fetch/CORS quirks. `index.html` loads `dsp.js` via a
  relative `<script src>`, so keep the two files together.
- **DSP unit tests**: `node tests/dsp.test.mjs` (zero deps, exits non-zero on failure). These cover
  the `dsp.js` primitives, including a guard that the fast `coherenceFromSegments` path stays
  numerically equal to the reference `welchCoherence`. Run them after touching any DSP code.
- Beyond that, "testing" is manual: open the page, pick a data source, press **Run Analysis**, and
  inspect the five tabs.
- The only browser dependency is **Chart.js 4.4.1**, loaded from a CDN `<script>` tag. Chart.js
  is used for the line-chart tabs (Time Series, FFT, window detail); everything else is hand-drawn
  on `<canvas>` 2D contexts.
- Data sources (toolbar **Source**) — two give *genuinely independent* per-station series (so
  coherence is measured, not modeled): **USGS** (`fetchUSGS` → `parseUSGSSeries`, default) needs no
  account but only serves US observatories (foreign `STATIONS` are dropped); **SuperMAG**
  (`fetchSuperMag` → `parseSuperMagSeries`) is global but needs a free user ID (the `logon` param,
  from the **SuperMAG ID** field). Both fetch by IAGA code — so a station's `id` must be its real
  IAGA code (e.g. College is `CMO`, not `COL`) — and both may be blocked by CORS from a static page
  (then synthetic fallback / use a proxy). **NOAA SWPC** (`fetchBaseSignal`, see `URLS`) returns one
  base signal from which stations are *derived* (modeled). **Synthetic** generates locally. On any
  fetch failure the app falls back to synthetic; the header provenance badge shows
  **REAL · source · cadence** vs **SYNTHETIC**, so the UI always produces output.
- Toolbar settings + active stations are persisted to `localStorage` and encoded in the URL hash
  (`saveSettings`/`loadSettings`, precedence URL > localStorage > defaults), so a configuration is
  shareable/reproducible. **PNG**/**CSV** export buttons dump the current coherogram.

## Signal-processing pipeline (the core)

All DSP is implemented from scratch in plain JS — no DSP library. The math is the heart of the
app; understand `runAnalysis()` before changing anything, as it orchestrates the whole pipeline:

1. **Acquire per-station signals.** `runAnalysis` dispatches on **Source**: SuperMAG fetches a real
   independent series per station (`fetchSuperMag` → `parseSuperMagSeries`); NOAA derives all
   stations from one base Bz series (`fetchBaseSignal`) via per-station `delay`/`ampScale`/noise (a
   *modeled* artifact, not independent data); synthetic uses `generateSynthetic`. Real sources derive
   the **true sampling rate `fs`** from timestamps (`medianCadenceSec`) rather than assuming 1 Hz;
   synthetic stays at `fs=1` (its Pc3/4/5 tones are real Hz at 1 Hz). The band is clamped to Nyquist
   (`fs/2`) — e.g. 1-minute data only resolves ≤~8.3 mHz (Pc5/low-Pc4). Stations that return no data
   are dropped and the rest are trimmed to a common length, so the effective `activeList` (used by
   every downstream consumer) may be smaller than the requested set. `STATIONS` is a global IAGA-code
   list; the globe view recenters on the active stations' centroid (`stationCentroid`).
2. **First-difference** each station signal (`firstDifference`) to emphasize the ULF band.
3. **Sliding window** across time (`sliceLen`/`stepSize` from the toolbar). For each window, compute
   each station's windowed FFT segments **once** (`welchSegments`), then form every pair's
   magnitude-squared coherence from that cache (`coherenceFromSegments`) — avoiding recomputing a
   station's FFT once per pair. Average the pairs, keep frequency bins ≤ `maxFreq`. Result is
   `coherogramData[timeIdx][freqIdx]`. (`welchCoherence` remains the reference impl, used for the
   full-signal matrix, the click-detail panel, and the equivalence test.)
4. **Significance threshold** (`computeSigThreshold`): analytic γ² threshold for p<0.05 given the
   number of Welch segments K. Used for the optional significance mask and per-pair ✓ flags.
5. **Event detection** (`detectEvents`): flags contiguous time windows whose band-mean coherence
   exceeds an adaptive threshold (max of the sig threshold, mean+1.5σ, and 0.35), classifies each
   into Pc3/Pc4/Pc5 by peak frequency.
6. Also computes a full mean-coherence **matrix** (`matrixMean`) and per-station **FFT PSD**
   (`stationFFT`) for the matrix/globe/FFT tabs.

Supporting primitives all live in **`dsp.js`**: `fft` (in-place iterative Cooley–Tukey, **requires
power-of-2 length** — use `nextPow2`), `infernoRGB` (approximate inferno colormap, drives every
heatmap and the globe arcs), `mulberry32` (seeded RNG for reproducible synthetic data), plus
`welchCoherence`/`welchSegments`/`coherenceFromSegments`, `firstDifference`, `computeSigThreshold`,
`medianCadenceSec`, and the per-source parsers `parseSuperMagSeries` / `parseUSGSSeries`
(records → `{values,times}`).
Keep `dsp.js` free of DOM/app state; if you add a primitive there, export it at the bottom and add a
test in `tests/dsp.test.mjs`.

## State & rendering model

- `lastResult` is the single global holding the full output of the most recent `runAnalysis()`.
  Every render function (`renderCoherogram`, `renderTimeSeries`, `renderFFTSpectra`, `renderMatrix`,
  `startGlobe`, `openDetail`) reads from it. After changing the analysis output shape, update every
  consumer.
- Tabs are switched with `switchTab(i)`; only the active tab is (re)rendered, and the globe's
  `requestAnimationFrame` loop is started/stopped on enter/leave to avoid background CPU use.
- Chart.js instances (`tsChart`, `fftChart`, `detailChart`) are destroyed and recreated on each
  render. Canvases are swapped via `freshCanvas(id)` to avoid Chart.js "canvas already in use"
  errors — keep that pattern when adding new charts.
- The coherogram, correlation matrix, colorbar, and globe are drawn imperatively on canvas (no
  Chart.js). The coherogram fills an `ImageData` pixel buffer directly for speed.
- `runAnalysis()` is `async` and yields to the UI via `await tick()` between heavy stages so the
  spinner/status updates render; preserve those yields when adding stages to long loops.
- **Auto-refresh** (`toggleAutoRefresh`/`scheduleRefresh`) re-runs the analysis every
  `REFRESH_INTERVAL` (60s) via chained `setTimeout`.

## Conventions

- Vanilla browser JS only, `'use strict'`, no modules/imports, no transpilation — keep it ES that
  runs directly in a modern browser. Typed arrays (`Float64Array`) are used throughout the DSP for
  performance; match that.
- Numeric inputs are clamped on read inside `runAnalysis` (e.g. `sliceLen` 120–3600); keep new
  inputs clamped the same way rather than trusting the DOM value.
- Section banners (`// ═══…`) organize the single script block; add new code under a matching
  banner rather than scattering it.
