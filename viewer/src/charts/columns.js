// How many columns a row of n charts gets when at most `most` fit side by side: rows of equal length, so that four
// charts stand four across or two by two (never three and one), six stand three by three, ten five by five.
// A plain function, checked without a browser (verify/units.mjs); charts/Grid.jsx states its answers to the style sheet.
export const columnsFor = (n, most) => { if (!(n > 0)) return 1; const rows = Math.ceil(n / Math.max(1, most)); return Math.ceil(n / rows); };
