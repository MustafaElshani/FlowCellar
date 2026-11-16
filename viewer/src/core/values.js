// How one stored number is written, by its id. Plain functions (no markup), so that they can be checked on their
// own (verify/units.mjs): the unit of a number must follow its id, whatever other words the id holds.
import * as F from './format.js';
import { FUNNEL_LABELS } from './derive.js';

export function fmtValue(id, v) {
  if (v == null) return F.DASH;
  if (typeof v === 'string') return id === 'funnel_usable_stage' && FUNNEL_LABELS[v] ? FUNNEL_LABELS[v].charAt(0).toLowerCase() + FUNNEL_LABELS[v].slice(1) : v.replace(/_/g, ' ').toLowerCase();
  // the unit an id ends in comes first: "fc_duty_half_life_h" is hours although it holds the word "duty"
  if (/_mv_per_h$/.test(id)) return v.toFixed(2) + ' mV per hour';
  if (/_per_h$/.test(id)) return F.sig(v) + ' Gb per hour';
  if (/_h$/.test(id)) return F.hours(v);
  if (/_c$/.test(id)) return v.toFixed(1) + ' °C';
  if (/_mv$/.test(id)) return Math.round(v) + ' mV';
  if (/_s$/.test(id)) return F.sig(v) + ' s';
  if (/_nt$|polya_median/.test(id)) return String(+v.toPrecision(4)) + ' nt';          // a poly(A) tail length, in nucleotides
  if (/^qual_/.test(id)) return F.q(v);
  if (/identity/.test(id)) return (v * 100).toFixed(1) + ' %';
  if (/_frac$|occupancy|duty/.test(id)) return F.pct(v);
  if (/gini/.test(id)) return v.toFixed(2);
  if (/speed/.test(id)) return F.sig(v) + ' bases/s';
  if (/_mb/.test(id)) return F.sig(v) + ' Mb';
  if (/_kb/.test(id)) return F.sig(v) + ' kb';
  if (/bases|n50|len_.*median$|mean$/.test(id)) return F.bases(v);
  if (/reads/.test(id)) return F.count(v);
  return Math.abs(v) >= 1000 ? F.int(v) : String(+v.toPrecision(4));
}
