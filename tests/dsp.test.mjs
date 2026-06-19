// Zero-dependency unit tests for the DSP primitives in ../dsp.js
// Run: node tests/dsp.test.mjs   (exits non-zero on any failure)
import { createRequire } from 'module';
const require = createRequire(import.meta.url);
const dsp = require('../dsp.js');
const { fft, nextPow2, welchCoherence, welchSegments, coherenceFromSegments,
        firstDifference, computeSigThreshold, mulberry32, medianCadenceSec,
        parseSuperMagSeries, parseUSGSSeries } = dsp;

let passed = 0, failed = 0;
function ok(cond, msg){ if(cond){ passed++; } else { failed++; console.error('  ✗ '+msg); } }
function approx(a, b, tol, msg){ ok(Math.abs(a-b) <= tol, `${msg} (got ${a}, want ${b}±${tol})`); }
function section(name){ console.log('• '+name); }

// ── naive DFT reference ──
function dft(re){
  const N=re.length, oRe=new Array(N).fill(0), oIm=new Array(N).fill(0);
  for(let k=0;k<N;k++) for(let n=0;n<N;n++){
    const a=-2*Math.PI*k*n/N; oRe[k]+=re[n]*Math.cos(a); oIm[k]+=re[n]*Math.sin(a);
  }
  return {oRe,oIm};
}

// ── nextPow2 ──
section('nextPow2');
ok(nextPow2(1)===1,'nextPow2(1)=1');
ok(nextPow2(3)===4,'nextPow2(3)=4');
ok(nextPow2(128)===128,'nextPow2(128)=128');
ok(nextPow2(129)===256,'nextPow2(129)=256');

// ── fft ──
section('fft');
{ // impulse → flat unit magnitude across all bins
  const N=16, re=new Float64Array(N), im=new Float64Array(N); re[0]=1;
  fft(re,im);
  let flat=true; for(let k=0;k<N;k++){ if(Math.abs(Math.hypot(re[k],im[k])-1)>1e-9) flat=false; }
  ok(flat,'impulse → flat magnitude 1');
}
{ // pure sinusoid at bin 3 → energy concentrated in bins 3 and N-3
  const N=32, b=3, re=new Float64Array(N), im=new Float64Array(N);
  for(let n=0;n<N;n++) re[n]=Math.cos(2*Math.PI*b*n/N);
  fft(re,im);
  let peak=0,peakK=-1,total=0;
  for(let k=0;k<N;k++){ const m=re[k]*re[k]+im[k]*im[k]; total+=m; if(m>peak){peak=m;peakK=k;} }
  ok(peakK===b||peakK===N-b,'sinusoid peak at expected bin');
  ok((peak/total)>0.45,'sinusoid energy concentrated at the peak');
}
{ // match naive DFT for a small random signal
  const N=8, re=new Float64Array(N), im=new Float64Array(N), rng=mulberry32(42);
  const x=Array.from({length:N},()=>rng()-0.5);
  for(let i=0;i<N;i++) re[i]=x[i];
  fft(re,im);
  const ref=dft(x);
  let maxErr=0; for(let k=0;k<N;k++) maxErr=Math.max(maxErr,Math.abs(re[k]-ref.oRe[k]),Math.abs(im[k]-ref.oIm[k]));
  ok(maxErr<1e-9,'fft matches naive DFT (maxErr '+maxErr.toExponential(2)+')');
}
{ // Parseval: sum|x|^2 == (1/N) sum|X|^2
  const N=16, re=new Float64Array(N), im=new Float64Array(N), rng=mulberry32(7);
  let timeE=0; for(let i=0;i<N;i++){ re[i]=rng()-0.5; timeE+=re[i]*re[i]; }
  fft(re,im);
  let freqE=0; for(let k=0;k<N;k++) freqE+=re[k]*re[k]+im[k]*im[k];
  approx(timeE, freqE/N, 1e-9, 'Parseval energy conservation');
}

// ── welchCoherence ──
section('welchCoherence');
{
  const nperseg=128, fs=4, n=1024, rng=mulberry32(99);
  const x=new Float64Array(n); for(let i=0;i<n;i++) x[i]=Math.sin(2*Math.PI*0.05*i)+0.3*(rng()-0.5);
  const r=welchCoherence(x,x,fs,nperseg);
  let allHigh=true; for(let k=1;k<r.coherence.length;k++) if(r.coherence[k]<0.999) allHigh=false;
  ok(allHigh,'identical signals → γ²≈1');
  approx(r.freqs[1], fs/nextPow2(nperseg), 1e-12, 'freqs[1] == fs/nfft');
}
{
  const nperseg=128, fs=1, n=4096, a=mulberry32(1), b=mulberry32(2);
  const x=new Float64Array(n), y=new Float64Array(n);
  for(let i=0;i<n;i++){ x[i]=a()-0.5; y[i]=b()-0.5; }
  const r=welchCoherence(x,y,fs,nperseg);
  let mean=0; for(let k=1;k<r.coherence.length;k++) mean+=r.coherence[k]; mean/=(r.coherence.length-1);
  ok(mean<0.35,'independent noise → low mean γ² (got '+mean.toFixed(3)+')');
}
ok(welchCoherence(new Float64Array(10),new Float64Array(10),1,128)===null,'too-short signal → null');

// ── welchSegments / coherenceFromSegments equivalence (perf path guard) ──
section('coherenceFromSegments ≈ welchCoherence');
{
  const nperseg=128, fs=2, n=900, ra=mulberry32(11), rb=mulberry32(22);
  const x=new Float64Array(n), y=new Float64Array(n);
  for(let i=0;i<n;i++){ x[i]=Math.sin(2*Math.PI*0.03*i)+0.5*(ra()-0.5); y[i]=Math.sin(2*Math.PI*0.03*i+0.4)+0.5*(rb()-0.5); }
  const ref=welchCoherence(x,y,fs,nperseg);
  const fast=coherenceFromSegments(welchSegments(x,fs,nperseg),welchSegments(y,fs,nperseg));
  let maxErr=0; for(let k=0;k<ref.coherence.length;k++) maxErr=Math.max(maxErr,Math.abs(ref.coherence[k]-fast.coherence[k]));
  ok(maxErr<1e-9,'optimized coherence matches reference (maxErr '+maxErr.toExponential(2)+')');
  ok(fast.nSeg===ref.nSeg,'same segment count');
}

// ── firstDifference ──
section('firstDifference');
{
  const d=firstDifference(new Float64Array([1,3,6,10]));
  ok(d.length===3 && d[0]===2 && d[1]===3 && d[2]===4,'firstDifference([1,3,6,10])=[2,3,4]');
}

// ── computeSigThreshold ──
section('computeSigThreshold');
{
  const t=computeSigThreshold(600,128);
  ok(t>0 && t<1,'threshold in (0,1)');
  ok(computeSigThreshold(600,128,0.01) > computeSigThreshold(600,128,0.05),'smaller α → higher threshold');
  ok(computeSigThreshold(3600,128) < computeSigThreshold(600,128),'more segments → lower threshold');
}

// ── mulberry32 determinism ──
section('mulberry32');
{
  const a=mulberry32(123), b=mulberry32(123);
  ok(a()===b() && a()===b(),'same seed → same stream');
  ok(mulberry32(1)()!==mulberry32(2)(),'different seeds → different output');
}

// ── medianCadenceSec (fs derivation for real data) ──
section('medianCadenceSec');
{
  // 1-minute-cadence timestamps, mag-7-day style "YYYY-MM-DD HH:MM:SS.sss"
  const tags=[]; const t0=Date.parse('2024-06-01 00:00:00.000');
  for(let i=0;i<10;i++) tags.push(t0 + i*60000);
  const cadence=medianCadenceSec(tags);
  approx(cadence,60,1e-6,'1-min tags → 60 s cadence');
  approx(1/cadence,1/60,1e-9,'→ fs ≈ 0.0167 Hz (Nyquist ≈ 8.3 mHz)');
  ok(medianCadenceSec([1000])===0,'single timestamp → 0 (caller falls back)');
  // robust to an out-of-order / duplicate glitch
  approx(medianCadenceSec([0,60000,60000,120000,180000]),60,1e-6,'median ignores a zero-gap glitch');
}

// ── parseSuperMagSeries (real per-station SuperMAG records) ──
section('parseSuperMagSeries');
{
  // documented shape: tval (unix s) + N/E/Z as {nez,geo}
  const recs=[
    {tval:1717200000, N:{nez:12.5,geo:11.1}, E:{nez:-3.1}, Z:{nez:4.0}},
    {tval:1717200060, N:{nez:13.0,geo:11.4}, E:{nez:-2.9}, Z:{nez:4.2}},
    {tval:1717200120, N:{nez:12.8}, E:{nez:-3.0}, Z:{nez:4.1}},
    {tval:1717200180, N:{nez:999999,geo:999999}, E:{nez:0}, Z:{nez:0}}, // missing → dropped
  ];
  const r=parseSuperMagSeries(recs,'N');
  ok(r.values.length===3,'drops the 999999 missing-data sentinel');
  ok(r.values[0]===12.5 && r.values[2]===12.8,'reads N.nez component');
  approx(medianCadenceSec(r.times),60,1e-6,'SuperMAG tval → 60 s cadence');
}
{
  // plain-number component form and a parseable date string
  const recs=[{tval:'2024-06-01T00:00:00',N:5},{tval:'2024-06-01T00:01:00',N:6}];
  const r=parseSuperMagSeries(recs,'N');
  ok(r.values.length===2 && r.values[1]===6,'handles numeric component + ISO tval');
}
ok(parseSuperMagSeries(null,'N').values.length===0,'non-array input → empty series');
ok(parseSuperMagSeries([{tval:1,N:{nez:7}}],'N').values[0]===7,'single record parsed');

// ── parseUSGSSeries (USGS Geomagnetism "Timeseries" JSON) ──
section('parseUSGSSeries');
{
  const json={
    type:'Timeseries',
    times:['2024-06-01T00:00:00Z','2024-06-01T00:01:00Z','2024-06-01T00:02:00Z','2024-06-01T00:03:00Z'],
    values:[
      {id:'H', values:[20512.3, 20512.5, 20512.1, null]}, // null gap dropped, large absolute kept
      {id:'Z', values:[48000, 48001, 48002, 48003]},
    ],
  };
  const r=parseUSGSSeries(json,'H');
  ok(r.values.length===3,'drops the null gap, keeps large absolute nT');
  ok(r.values[0]===20512.3 && r.values[2]===20512.1,'selects the requested element (H)');
  approx(medianCadenceSec(r.times),60,1e-6,'USGS ISO times → 60 s cadence');
  ok(parseUSGSSeries(json,'Z').values[1]===48001,'can select a different element');
}
ok(parseUSGSSeries(null,'H').values.length===0,'non-Timeseries input → empty');
ok(parseUSGSSeries({times:['2024-06-01T00:00:00Z'],values:[{id:'X',values:[1]}]},'H').values[0]===1,
  'falls back to first series when element id is absent');

console.log(`\n${passed} passed, ${failed} failed`);
process.exit(failed?1:0);
