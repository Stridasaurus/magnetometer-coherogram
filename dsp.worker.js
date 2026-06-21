'use strict';
importScripts('./dsp.js');

self.addEventListener('message', e => {
  if (e.data.type === 'analyze') runDSP(e.data);
  // 'cancel' is handled by worker.terminate() on the main thread — no flag needed
});

function runDSP({stationData, settings}) {
  const {activeList, fs, cadenceSec, sliceSamples, stepSamples, maxFreqHz,
         npersegBase, skipPairCoh, useSynth, src, pickDate,
         fillFrac, warnStations, requestedCount} = settings;

  // First-difference each station signal
  const diffSigs = {};
  for (const s of activeList) diffSigs[s.id] = firstDifference(stationData[s.id]);
  const nSamples = diffSigs[activeList[0].id].length;
  const nPairs = activeList.length * (activeList.length - 1) / 2;

  // Clamp nperseg so at least 2 Welch segments fit in sliceSamples
  let nperseg = npersegBase;
  while (nperseg > 4 && nperseg * 1.5 > sliceSamples) nperseg >>= 1;
  const sigThreshold = computeSigThreshold(sliceSamples, nperseg);

  // Build frequency grid from a test coherence call
  const testRes = welchCoherence(
    diffSigs[activeList[0].id].slice(0, sliceSamples),
    diffSigs[activeList[1].id].slice(0, sliceSamples),
    fs, nperseg);
  if (!testRes) { self.postMessage({type:'error', msg:'Signal too short for coherence'}); return; }
  const allFreqs = testRes.freqs;
  const ulfMask = [];
  // Start at k=1 to exclude the DC bin (0 Hz) — unstable after first-differencing
  for (let k = 1; k < allFreqs.length; k++) if (allFreqs[k] <= maxFreqHz) ulfMask.push(k);
  const ulfFreqs = ulfMask.map(k => allFreqs[k]);
  const nULF = ulfFreqs.length;
  const totalSlices = Math.ceil((nSamples - sliceSamples) / stepSamples + 1);

  // Sliding window coherogram
  // pairMatSum accumulates windowed coherence per pair for matrix computation —
  // gives a consistent "temporal mean of windowed γ²" at all network sizes.
  const times = [], windowStarts = [], coherogramData = [], coherogramPhase = [], pairCohMap = {};
  const pairMatSum = {};
  let sliceIdx = 0;
  for (let start = 0; start + sliceSamples <= nSamples; start += stepSamples) {
    times.push((start + sliceSamples / 2) / fs);  // seconds, not sample indices
    windowStarts.push(start);
    // Compute each station's FFT segments once, reuse across all pairs
    const segs = {};
    for (const s of activeList) segs[s.id] = welchSegments(diffSigs[s.id].slice(start, start + sliceSamples), fs, nperseg);
    const sumCoh = new Float64Array(nULF);
    const sumPhSin = new Float64Array(nULF), sumPhCos = new Float64Array(nULF);
    let pCnt = 0;
    for (let i = 0; i < activeList.length; i++) for (let j = i + 1; j < activeList.length; j++) {
      const ai = activeList[i], bi = activeList[j];
      const res = coherenceFromSegments(segs[ai.id], segs[bi.id]);
      if (res) {
        const pk = `${ai.id}:${bi.id}`, pRow = new Float64Array(nULF);
        for (let u = 0; u < nULF; u++) {
          const c = res.coherence[ulfMask[u]];
          sumCoh[u] += c;
          sumPhSin[u] += Math.sin(res.phase[ulfMask[u]]);
          sumPhCos[u] += Math.cos(res.phase[ulfMask[u]]);
          pRow[u] = c;
        }
        // Accumulate for matrix (always — gives consistent definition regardless of skipPairCoh)
        if (!pairMatSum[pk]) pairMatSum[pk] = {sum: new Float64Array(nULF), cnt: 0};
        for (let u = 0; u < nULF; u++) pairMatSum[pk].sum[u] += pRow[u];
        pairMatSum[pk].cnt++;
        // Store full per-window rows only for small networks (per-pair coherogram modal)
        if (!skipPairCoh) { if (!pairCohMap[pk]) pairCohMap[pk] = []; pairCohMap[pk].push(pRow); }
        pCnt++;
      }
    }
    const row = new Float64Array(nULF), phRow = new Float64Array(nULF);
    if (pCnt) for (let u = 0; u < nULF; u++) {
      row[u] = sumCoh[u] / pCnt;
      phRow[u] = Math.atan2(sumPhSin[u], sumPhCos[u]);
    }
    coherogramData.push(row);
    coherogramPhase.push(phRow);
    sliceIdx++;
    if (sliceIdx % 8 === 0) self.postMessage({type:'progress', msg:`Slice ${sliceIdx}/${totalSlices}…`});
  }
  const nTimes = times.length;

  // Cross-coherence matrix — temporal mean of windowed γ² via pairMatSum
  // (same definition for all network sizes, no welchCoherence fallback needed)
  const matrixMean = activeList.map(() => new Float64Array(activeList.length));
  for (let i = 0; i < activeList.length; i++) for (let j = 0; j < activeList.length; j++) {
    if (i === j) { matrixMean[i][j] = 1; continue; }
    const [a, b] = i < j ? [activeList[i], activeList[j]] : [activeList[j], activeList[i]];
    const pk = `${a.id}:${b.id}`, entry = pairMatSum[pk];
    if (entry && entry.cnt) {
      let bm = 0; for (let u = 0; u < nULF; u++) bm += entry.sum[u] / nULF;
      matrixMean[i][j] = bm / entry.cnt;
    }
  }

  // Per-station FFT power spectral density
  const stationFFT = {};
  for (const s of activeList) {
    const sig = diffSigs[s.id];
    const nfft = nextPow2(Math.min(sig.length, 4096));
    const re = new Float64Array(nfft), im = new Float64Array(nfft);
    for (let i = 0; i < nfft && i < sig.length; i++) re[i] = sig[i] * 0.5 * (1 - Math.cos(2 * Math.PI * i / (nfft - 1)));
    fft(re, im);
    const half = (nfft >> 1) + 1, psd = new Float64Array(half), pFreqs = new Float64Array(half);
    for (let k = 0; k < half; k++) { psd[k] = Math.log10(Math.max(1e-30, re[k]*re[k]+im[k]*im[k])); pFreqs[k] = k * fs / nfft; }
    stationFFT[s.id] = {psd, freqs: pFreqs};
  }

  // Peak coherence stats
  let maxCoh = 0, peakFreqHz = 0;
  for (let t = 0; t < nTimes; t++) for (let u = 0; u < nULF; u++) {
    if (coherogramData[t][u] > maxCoh) { maxCoh = coherogramData[t][u]; peakFreqHz = ulfFreqs[u]; }
  }

  self.postMessage({type:'result', result:{
    coherogramData, coherogramPhase, pairCohMap, times, windowStarts, ulfFreqs, ulfMask, nTimes, nULF,
    activeList, diffSigs, stationFFT, matrixMean, maxCoh, peakFreqHz, nSamples, nPairs,
    useSynth, src, sliceSamples, sigThreshold, allFreqs, nperseg, fs, cadenceSec,
    pickDate, fillFrac, warnStations, requestedCount
  }});
}
