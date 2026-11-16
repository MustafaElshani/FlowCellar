// Where each channel sits on the flow cell map.
//
// A bundle names only the layout ("promethion_3000", "minion_512"); no layout table from the manufacturer ships
// with FlowCellar, so no map here claims true positions:
//
// promethion_3000: 12 blocks of 250 channels side by side, each block 10 columns x 25 rows, channels running left
//   to right and then down inside a block: 120 columns x 25 rows. This arrangement is INFERRED from how channels
//   behave in real runs (channels that are neighbours under this rule behave alike in both directions, which no
//   other simple arrangement of the 3,000 numbers gives); it is not taken from the manufacturer. Which corner
//   channel 1 sits in, and whether blocks are mirrored, is not known. The map says "positions inferred".
//
// minion_512: not known. Channels are drawn in number order, 32 per row, and the map says so.
export const LAYOUTS = {
  promethion_3000: {
    cols: 120, rows: 25, physical: 'inferred',
    note: 'Channel positions are inferred from how neighbouring channels behave in real runs; they are not taken from the manufacturer’s layout table, so corners and mirroring may differ from the flow cell.',
    pos: ch => { const b = (ch / 250) | 0, i = ch % 250; return [(i / 10) | 0, b * 10 + i % 10]; },
    at: (r, c) => { const b = (c / 10) | 0; return b * 250 + r * 10 + c % 10; },
  },
  minion_512: {
    cols: 32, rows: 16, physical: 'order',
    note: 'Channels are drawn in number order, 32 per row: the physical layout of this flow cell type is not known to FlowCellar, so neighbours here are not neighbours on the flow cell.',
    pos: ch => [(ch / 32) | 0, ch % 32],
    at: (r, c) => r * 32 + c,
  },
};
export function layoutFor(name, n) {
  if (LAYOUTS[name]) return LAYOUTS[name];
  const cols = Math.ceil(Math.sqrt(n * 2)), rows = Math.ceil(n / cols);
  return { cols, rows, physical: 'order', note: 'Channels are drawn in number order.', pos: ch => [(ch / cols) | 0, ch % cols], at: (r, c) => { const ch = r * cols + c; return ch < n ? ch : null; } };
}
