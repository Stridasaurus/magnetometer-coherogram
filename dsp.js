'use strict';
// ═══════════════════════════════════════════════════════════
//  DSP PRIMITIVES
//  Pure, dependency-free signal-processing functions shared by
//  index.html (loaded as a plain browser <script>, functions
//  become globals) and tests/dsp.test.mjs (Node, via require).
//  Keep this file free of DOM / app state.
// ═══════════════════════════════════════════════════════════

// ── Inferno colormap (approximate) ──
function infernoRGB(t) {
  t = Math.max(0, Math.min(1, t));
  const r = 255*(0.001 + t*(0.8*(1-Math.exp(-5*t)) + 0.8*t*t));
  const g = t < 0.5 ? 255*(0.05+t*0.4) : 255*(0.25+(t-0.5)*2*0.65);
  const b = t < 0.3 ? 255*(0.5*t/0.3) : t < 0.7 ? 255*(0.5*(1-(t-0.3)/0.4)) : 255*0.02;
  return [Math.max(0,Math.min(255,Math.round(r))), Math.max(0,Math.min(255,Math.round(g))), Math.max(0,Math.min(255,Math.round(b)))];
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
  const coherence=new Float64Array(half), freqs=new Float64Array(half);
  for (let k=0;k<half;k++){
    freqs[k]=k*fs/nfft;
    const num=PxyRe[k]*PxyRe[k]+PxyIm[k]*PxyIm[k];
    const den=Pxx[k]*Pyy[k];
    coherence[k]=den>1e-30?Math.min(1,num/den):0;
  }
  return {freqs, coherence, nSeg};
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
  const coherence=new Float64Array(half);
  for (let k=0;k<half;k++){
    const num=PxyRe[k]*PxyRe[k]+PxyIm[k]*PxyIm[k];
    const den=a.Pxx[k]*b.Pxx[k];
    coherence[k]=den>1e-30?Math.min(1,num/den):0;
  }
  return {freqs:a.freqs, coherence, nSeg};
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

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { infernoRGB, fft, nextPow2, welchCoherence, welchSegments,
    coherenceFromSegments, firstDifference, computeSigThreshold, mulberry32, medianCadenceSec };
}
