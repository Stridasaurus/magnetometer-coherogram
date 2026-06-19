'use strict';
// ═══════════════════════════════════════════════════════════
//  DSP PRIMITIVES
//  Pure, dependency-free signal-processing functions shared by
//  index.html (loaded as a plain browser <script>, functions
//  become globals) and tests/dsp.test.mjs (Node, via require).
//  Keep this file free of DOM / app state.
// ═══════════════════════════════════════════════════════════

// ── Inferno colormap — 256-entry LUT from matplotlib's official data ──
const _IR=new Uint8Array([0,1,1,1,2,2,2,3,4,4,5,6,7,8,9,10,11,12,13,14,16,17,18,20,21,22,24,25,27,28,30,31,33,35,36,38,40,41,43,45,47,49,50,52,54,56,57,59,61,62,64,66,68,69,71,73,74,76,77,79,81,82,84,85,87,89,90,92,93,95,97,98,100,101,103,105,106,108,109,111,113,114,116,117,119,120,122,124,125,127,128,130,132,133,135,136,138,140,141,143,144,146,147,149,151,152,154,155,157,159,160,162,163,165,166,168,169,171,173,174,176,177,179,180,182,183,185,186,188,189,191,192,193,195,196,198,199,200,202,203,204,206,207,208,210,211,212,213,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,235,236,237,238,239,239,240,241,241,242,243,243,244,245,245,246,246,247,247,248,248,248,249,249,249,250,250,250,251,251,251,251,251,252,252,252,252,252,252,252,252,252,252,252,252,251,251,251,251,251,250,250,250,250,249,249,249,248,248,247,247,246,246,245,245,244,244,244,243,243,242,242,242,241,241,241,241,242,242,243,243,244,245,246,248,249,250,252]);
const _IG=new Uint8Array([0,0,1,1,1,2,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10,11,11,11,12,12,12,12,12,12,12,12,12,12,11,11,11,11,10,10,10,10,9,9,9,9,9,9,10,10,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,18,18,19,19,20,21,21,22,22,23,24,24,25,25,26,26,27,28,28,29,29,30,30,31,32,32,33,33,34,34,35,35,36,37,37,38,38,39,39,40,41,41,42,42,43,44,44,45,46,46,47,48,48,49,50,50,51,52,53,53,54,55,56,57,58,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,74,75,76,77,78,80,81,82,83,85,86,87,89,90,92,93,94,96,97,99,100,102,103,105,106,108,110,111,113,115,116,118,120,121,123,125,126,128,130,132,133,135,137,139,140,142,144,146,148,150,151,153,155,157,159,161,163,165,166,168,170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,199,201,203,205,207,209,211,213,215,217,219,221,223,225,227,229,230,232,234,236,237,239,241,242,244,245,246,248,249,250,251,252,253,255]);
const _IB=new Uint8Array([4,5,6,8,10,12,14,16,18,20,23,25,27,29,31,34,36,38,41,43,45,48,50,52,55,57,60,62,65,67,69,72,74,76,79,81,83,85,87,89,91,92,94,95,97,98,99,100,101,102,103,104,104,105,106,106,107,107,108,108,108,109,109,109,110,110,110,110,110,110,110,110,110,110,110,110,110,110,110,110,110,110,110,110,109,109,109,109,109,108,108,108,107,107,107,106,106,105,105,105,104,104,103,103,102,102,101,100,100,99,99,98,97,96,96,95,94,94,93,92,91,90,90,89,88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,70,69,68,67,66,65,63,62,61,60,59,58,56,55,54,53,52,51,49,48,47,46,45,43,42,41,40,38,37,36,35,33,32,31,29,28,27,25,24,23,21,20,19,18,16,15,14,12,11,10,9,8,7,7,6,6,6,6,7,7,8,9,10,12,13,15,17,18,20,22,24,26,29,31,33,35,38,40,42,45,47,50,53,55,58,61,64,67,70,73,76,79,83,86,90,93,97,101,105,109,113,117,121,125,130,134,138,142,146,150,154,157,161,164]);
function infernoRGB(t) {
  t = Math.max(0, Math.min(1, t));
  const s = t*255, i = s|0, j = Math.min(255,i+1), f = s-i;
  return [Math.round(_IR[i]+(_IR[j]-_IR[i])*f), Math.round(_IG[i]+(_IG[j]-_IG[i])*f), Math.round(_IB[i]+(_IB[j]-_IB[i])*f)];
}

// ── In-place iterative Cooley–Tukey FFT (requires power-of-2 length) ──
function fft(re, im) {
  const N = re.length;
  let j = 0;
  for (let i = 1; i < N; i++) {
    let bit = N >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { [re[i],re[j]]=[re[j],re[i]]; [im[i],im[j]]=[im[j],im[i]]; }
  }
  for (let len = 2; len <= N; len <<= 1) {
    const ang = -2*Math.PI/len, wRe = Math.cos(ang), wIm = Math.sin(ang);
    for (let i = 0; i < N; i += len) {
      let cRe = 1, cIm = 0;
      for (let k = 0; k < len/2; k++) {
        const uRe=re[i+k], uIm=im[i+k];
        const vRe=re[i+k+len/2]*cRe - im[i+k+len/2]*cIm;
        const vIm=re[i+k+len/2]*cIm + im[i+k+len/2]*cRe;
        re[i+k]=uRe+vRe; im[i+k]=uIm+vIm;
        re[i+k+len/2]=uRe-vRe; im[i+k+len/2]=uIm-vIm;
        const tmp=cRe*wRe-cIm*wIm; cIm=cRe*wIm+cIm*wRe; cRe=tmp;
      }
    }
  }
}
function nextPow2(n) { let p=1; while(p<n) p<<=1; return p; }

// ── Welch magnitude-squared coherence (reference implementation) ──
function welchCoherence(x, y, fs, nperseg) {
  const N = nperseg, step = N>>1, nfft = nextPow2(N), half = (nfft>>1)+1;
  const win = new Float64Array(N);
  for (let i=0;i<N;i++) win[i]=0.5*(1-Math.cos(2*Math.PI*i/(N-1)));
  const Pxx=new Float64Array(half), Pyy=new Float64Array(half);
  const PxyRe=new Float64Array(half), PxyIm=new Float64Array(half);
  let nSeg=0;
  for (let s=0; s+N<=x.length; s+=step) {
    const rX=new Float64Array(nfft),iX=new Float64Array(nfft);
    const rY=new Float64Array(nfft),iY=new Float64Array(nfft);
    for (let i=0;i<N;i++){rX[i]=x[s+i]*win[i]; rY[i]=y[s+i]*win[i];}
    fft(rX,iX); fft(rY,iY);
    for (let k=0;k<half;k++){
      Pxx[k]+=rX[k]*rX[k]+iX[k]*iX[k];
      Pyy[k]+=rY[k]*rY[k]+iY[k]*iY[k];
      PxyRe[k]+=rX[k]*rY[k]+iX[k]*iY[k];
      PxyIm[k]+=iX[k]*rY[k]-rX[k]*iY[k];
    }
    nSeg++;
  }
  if (!nSeg) return null;
  const coherence=new Float64Array(half), freqs=new Float64Array(half), phase=new Float64Array(half);
  for (let k=0;k<half;k++){
    freqs[k]=k*fs/nfft;
    const num=PxyRe[k]*PxyRe[k]+PxyIm[k]*PxyIm[k];
    const den=Pxx[k]*Pyy[k];
    coherence[k]=den>1e-30?Math.min(1,num/den):0;
    phase[k]=Math.atan2(PxyIm[k],PxyRe[k]);
  }
  return {freqs, coherence, phase, nSeg};
}

// ── Per-station Welch segment spectra (cache for fast pairwise coherence) ──
// Computes each Hann-windowed segment's complex spectrum once for a single
// signal, plus the accumulated auto-spectrum (Pxx). Uses the *exact* window,
// nfft and stepping as welchCoherence so coherenceFromSegments reproduces it.
function welchSegments(x, fs, nperseg) {
  const N = nperseg, step = N>>1, nfft = nextPow2(N), half = (nfft>>1)+1;
  const win = new Float64Array(N);
  for (let i=0;i<N;i++) win[i]=0.5*(1-Math.cos(2*Math.PI*i/(N-1)));
  const segRe=[], segIm=[];
  const Pxx=new Float64Array(half), freqs=new Float64Array(half);
  for (let k=0;k<half;k++) freqs[k]=k*fs/nfft;
  let nSeg=0;
  for (let s=0; s+N<=x.length; s+=step) {
    const rX=new Float64Array(nfft), iX=new Float64Array(nfft);
    for (let i=0;i<N;i++) rX[i]=x[s+i]*win[i];
    fft(rX,iX);
    const re=new Float64Array(half), im=new Float64Array(half);
    for (let k=0;k<half;k++){ re[k]=rX[k]; im[k]=iX[k]; Pxx[k]+=rX[k]*rX[k]+iX[k]*iX[k]; }
    segRe.push(re); segIm.push(im);
    nSeg++;
  }
  return {half, nSeg, freqs, segRe, segIm, Pxx};
}

// ── Magnitude-squared coherence from two cached welchSegments results ──
// Numerically equivalent to welchCoherence(a_signal, b_signal, ...).
function coherenceFromSegments(a, b) {
  const half=a.half, nSeg=a.nSeg;
  if (!nSeg) return null;
  const PxyRe=new Float64Array(half), PxyIm=new Float64Array(half);
  for (let s=0;s<nSeg;s++){
    const aRe=a.segRe[s], aIm=a.segIm[s], bRe=b.segRe[s], bIm=b.segIm[s];
    for (let k=0;k<half;k++){
      PxyRe[k]+=aRe[k]*bRe[k]+aIm[k]*bIm[k];
      PxyIm[k]+=aIm[k]*bRe[k]-aRe[k]*bIm[k];
    }
  }
  const coherence=new Float64Array(half), phase=new Float64Array(half);
  for (let k=0;k<half;k++){
    const num=PxyRe[k]*PxyRe[k]+PxyIm[k]*PxyIm[k];
    const den=a.Pxx[k]*b.Pxx[k];
    coherence[k]=den>1e-30?Math.min(1,num/den):0;
    phase[k]=Math.atan2(PxyIm[k],PxyRe[k]);
  }
  return {freqs:a.freqs, coherence, phase, nSeg};
}

// ── First difference (emphasise the ULF band) ──
function firstDifference(sig) {
  const out=new Float64Array(sig.length-1);
  for(let i=0;i<out.length;i++) out[i]=sig[i+1]-sig[i];
  return out;
}

// ── Analytic γ² significance threshold for p<alpha given K Welch segments ──
function computeSigThreshold(sliceLen, nperseg, alpha=0.05) {
  const K=Math.max(2, Math.floor((sliceLen-nperseg)/(nperseg>>1))+1);
  return 1-Math.pow(alpha, 1/(K-1));
}

// ── Seeded RNG (reproducible synthetic data) ──
function mulberry32(a){ return()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;}; }

// ── Median sample cadence (seconds) from an array of epoch-ms timestamps ──
// Used to derive the true fs of fetched data instead of assuming 1 Hz.
function medianCadenceSec(times){
  const d=[];
  for(let i=1;i<times.length;i++){const dt=times[i]-times[i-1]; if(isFinite(dt)&&dt>0) d.push(dt);}
  if(!d.length) return 0;
  d.sort((a,b)=>a-b);
  return d[Math.floor(d.length/2)]/1000;
}

// ── Parse one SuperMAG station's JSON records into a {values,times} series ──
// Defensive against the documented shapes: a record carries a timestamp
// (`tval` unix-seconds, or a parseable date string) and N/E/Z components that
// may be plain numbers or {nez,geo} objects. SuperMAG's ~999999 missing-data
// sentinel and any non-finite samples are dropped. comp defaults to 'N'.
function parseSuperMagSeries(records, comp){
  comp = comp || 'N';
  const values=[], times=[];
  if(!Array.isArray(records)) return {values, times};
  for(const rec of records){
    if(!rec || typeof rec!=='object') continue;
    let tms;
    if(typeof rec.tval==='number') tms=rec.tval*1000;
    else if(rec.tval!=null) tms=Date.parse(rec.tval);
    else if(rec.time_tag!=null) tms=Date.parse(rec.time_tag);
    else tms=NaN;
    let v=rec[comp];
    if(v && typeof v==='object') v=(v.nez!=null?v.nez:v.geo);
    v=parseFloat(v);
    if(isFinite(v) && Math.abs(v)<1e5 && isFinite(tms)){ values.push(v); times.push(tms); }
  }
  return {values, times};
}

// ── Parse the USGS Geomagnetism web-service "Timeseries" JSON ──
// Shape: { times:[ISO,…], values:[ { id:'H', values:[num|null,…] }, … ] }.
// Pairs the chosen element's samples with `times`, dropping null/non-finite
// gaps. Values are absolute field (nT, ~tens of thousands) — no magnitude
// filter; the first-difference downstream removes the DC level.
function parseUSGSSeries(json, element){
  const out={values:[], times:[]};
  if(!json || typeof json!=='object' || !Array.isArray(json.times) || !Array.isArray(json.values)) return out;
  const series=json.values.find(s=>s && s.id===element) || json.values[0];
  if(!series || !Array.isArray(series.values)) return out;
  const vals=series.values, times=json.times, n=Math.min(vals.length, times.length);
  for(let i=0;i<n;i++){
    const v=parseFloat(vals[i]), t=Date.parse(times[i]);
    if(isFinite(v) && isFinite(t)){ out.values.push(v); out.times.push(t); }
  }
  return out;
}

// ── Circular mean of an array of phase values (radians) ──
function circMeanPhase(phases) {
  let x=0, y=0;
  for(const p of phases){ x+=Math.cos(p); y+=Math.sin(p); }
  return Math.atan2(y, x);
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { infernoRGB, fft, nextPow2, welchCoherence, welchSegments,
    coherenceFromSegments, firstDifference, computeSigThreshold, mulberry32,
    medianCadenceSec, parseSuperMagSeries, parseUSGSSeries, circMeanPhase };
}
