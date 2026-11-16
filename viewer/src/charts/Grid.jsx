// A row of charts. Up to 1,680 px it is the two-column grid of the run page (three columns with wide = 3). On a wider
// window the row has as many columns as fit at 420 px, and the charts are spread over rows of equal length: four charts
// stand four across or two by two, never three and one; a chart left over in the last row keeps the width of the
// others. A column is at most --fc-chart-wide (styles/app.css), and every row starts at the left gutter.
//
// The number of charts is counted in the page itself (a child may render nothing), before the first paint; which of
// the stated column counts applies is decided by the style sheet from the width of the window, so nothing is measured.
import { useLayoutEffect, useRef } from 'preact/hooks';
import { columnsFor } from './columns.js';

export function Grid({ children, wide, cls }) {
  const ref = useRef(null);
  useLayoutEffect(() => {
    const el = ref.current; if (!el) return; const n = el.children.length;
    if (el.dataset.n === String(n)) return; el.dataset.n = n;
    for (let most = 1; most <= 8; most++) el.style.setProperty('--fc-k' + most, columnsFor(n, most));
  });
  return <div ref={ref} class={'fc-grid fc-grid--' + (wide || 2) + ' fc-grid--n' + (cls ? ' ' + cls : '')}>{children}</div>;
}
