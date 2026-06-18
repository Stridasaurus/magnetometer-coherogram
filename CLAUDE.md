# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-page, zero-build web app that visualizes a **network coherogram** of ground/space
magnetometer data to detect ULF (Pc3–Pc5, ~2–50 mHz) waves. The entire app — HTML, CSS, and
JavaScript — lives in **`index.html`**. There is no package.json, build step, bundler, test
runner, or backend.

## Running / developing

- Open `index.html` directly in a browser, or serve it (e.g. `python3 -m http.server`) and load
  it. Serving over HTTP avoids any `file://` fetch/CORS quirks.
- No install, build, lint, or test commands exist. "Testing" is manual: open the page, pick a
  data source, press **Run Analysis**, and inspect the five tabs.
- The only external dependency is **Chart.js 4.4.1**, loaded from a CDN `<script>` tag. Chart.js
  is used for the line-chart tabs (Time Series, FFT, window detail); everything else is hand-drawn
  on `<canvas>` 2D contexts.
- Live data is fetched from NOAA SWPC endpoints (see `URLS` in the script). These require network
  access; on any fetch failure the app silently falls back to a synthetic substorm signal, so the
  UI always produces output.

## Signal-processing pipeline (the core)

All DSP is implemented from scratch in plain JS — no DSP library. The math is the heart of the
app; understand `runAnalysis()` before changing anything, as it orchestrates the whole pipeline:

1. **Fetch** a single base Bz time series (`fetchBaseSignal`), or generate synthetic per-station
   signals (`generateSynthetic`). For real data, each of the 6 `STATIONS` is derived from the one
   base signal by applying a per-station time `delay`, `ampScale`, and added noise — so inter-station
   coherence is a modeled artifact, not independently measured data.
2. **First-difference** each station signal (`firstDifference`) to emphasize the ULF band.
3. **Sliding window** across time (`sliceLen`/`stepSize` from the toolbar). For each window, compute
   **Welch magnitude-squared coherence** (`welchCoherence`) for every active station pair, average
   the pairs, and keep only frequency bins ≤ `maxFreq`. The result is `coherogramData[timeIdx][freqIdx]`.
4. **Significance threshold** (`computeSigThreshold`): analytic γ² threshold for p<0.05 given the
   number of Welch segments K. Used for the optional significance mask and per-pair ✓ flags.
5. **Event detection** (`detectEvents`): flags contiguous time windows whose band-mean coherence
   exceeds an adaptive threshold (max of the sig threshold, mean+1.5σ, and 0.35), classifies each
   into Pc3/Pc4/Pc5 by peak frequency.
6. Also computes a full mean-coherence **matrix** (`matrixMean`) and per-station **FFT PSD**
   (`stationFFT`) for the matrix/globe/FFT tabs.

Supporting primitives: `fft` (in-place iterative Cooley–Tukey, **requires power-of-2 length** — use
`nextPow2`), `infernoRGB` (approximate inferno colormap, drives every heatmap and the globe arcs),
`mulberry32` (seeded RNG for reproducible synthetic data).

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
