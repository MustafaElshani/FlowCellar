// Theme: every colour a chart uses is read from the design tokens. No colour literal in chart code.
let cs = null, memo = new Map(), probe = null, gen = 0;
export function resetTheme() { cs = null; memo = new Map(); lut = null; gen++; }
// counts the changes of theme: a chart that was painted under an earlier count holds the colours of that theme
export const themeGen = () => gen;
export function tk(name) {
  let v = memo.get(name);
  if (v !== undefined) return v;
  if (!cs) cs = getComputedStyle(document.documentElement);
  v = cs.getPropertyValue(name).trim();
  memo.set(name, v);
  return v;
}
export const font = (weight) => (weight ? weight + ' ' : '') + (tk('--fc-chart-font') || '11px') + ' ' + (tk('--fc-font-sans') || 'sans-serif');

// any CSS colour -> [r, g, b] through a canvas (also resolves system colours under forced colours)
export function rgb(color) {
  const key = 'rgb:' + color; let v = memo.get(key); if (v) return v;
  if (!probe) { probe = document.createElement('canvas'); probe.width = probe.height = 1; }
  const g = probe.getContext('2d', { willReadFrequently: true });
  g.clearRect(0, 0, 1, 1); g.fillStyle = '#000'; g.fillStyle = color; g.fillRect(0, 0, 1, 1);
  const d = g.getImageData(0, 0, 1, 1).data; v = [d[0], d[1], d[2]]; memo.set(key, v); return v;
}
export function luminance(color) { const [r, g, b] = rgb(color).map(c => { c /= 255; return c <= 0.03928 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4); }); return 0.2126 * r + 0.7152 * g + 0.0722 * b; }
export const onData = color => luminance(color) > 0.36 ? tk('--fc-on-data-dark') : tk('--fc-on-data-light');
export function alpha(color, a) { const [r, g, b] = rgb(color); return `rgba(${r},${g},${b},${a})`; }

// ---- sequential ramp: the nine steps, interpolated in OKLab ----
const lin = c => { c /= 255; return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4); };
const gam = c => Math.round(255 * Math.max(0, Math.min(1, c <= 0.0031308 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055)));
function toLab([r, g, b]) {
  r = lin(r); g = lin(g); b = lin(b);
  const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b), m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b), s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
  return [0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s, 1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s, 0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s];
}
function fromLab([L, a, b]) {
  const l = Math.pow(L + 0.3963377774 * a + 0.2158037573 * b, 3), m = Math.pow(L - 0.1055613458 * a - 0.0638541728 * b, 3), s = Math.pow(L - 0.0894841775 * a - 1.291485548 * b, 3);
  return [gam(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s), gam(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s), gam(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s)];
}
let lut = null;
// 256 colours from seq-1 to seq-9; entry t in 0..1
export function seqLUT() {
  if (lut) return lut;
  const stops = []; for (let i = 1; i <= 9; i++) stops.push(toLab(rgb(tk('--fc-seq-' + i))));
  lut = [];
  for (let i = 0; i < 256; i++) {
    const t = i / 255 * 8, k = Math.min(7, Math.floor(t)), f = t - k, a = stops[k], b = stops[k + 1];
    lut.push(fromLab([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]));
  }
  return lut;
}
export const seqColor = t => { const c = seqLUT()[Math.max(0, Math.min(255, Math.round(t * 255)))]; return `rgb(${c[0]},${c[1]},${c[2]})`; };
