# Third-party notices

FlowCellar itself is under the licence in `LICENSE`. This file lists the work of others that is **inside**
the files FlowCellar ships, and the terms it comes under.

## What is bundled

The viewer (`flowcellar/viewer_dist/viewer.core.js`, and therefore every hub page and every report FlowCellar
writes) is built from FlowCellar's own code and from these three packages. Their exact versions are fixed in
`viewer/package-lock.json`.

| Package | Version | Licence | Copyright |
|---|---|---|---|
| [preact](https://github.com/preactjs/preact) | 11.2.0 | MIT | Copyright (c) 2015-present Jason Miller |
| [@preact/signals](https://github.com/preactjs/signals) | 2.11.4 | MIT | Copyright (c) 2022-present Preact Team |
| [@preact/signals-core](https://github.com/preactjs/signals) | 1.14.5 | MIT | Copyright (c) 2022-present Preact Team |

Nothing else is bundled. No font, icon set, chart library or stylesheet of another project is inside the
package: text uses the reader's system fonts, the charts are drawn by FlowCellar's own canvas code, and the
schema, the catalogue and the styles were written for this project.

The word "FlowCellar" in the logo was set in Inter Display (Copyright (c) The Inter Project Authors, SIL Open
Font License 1.1) and is delivered as outlines. No font file is shipped, and the Open Font License does not
apply to documents or artwork made with a font.

## What is not bundled

These are installed beside FlowCellar by the installer, under their own licences, and no part of them is
inside FlowCellar's files: numpy, pyarrow, click (always), pysam (extra `bam`), jsonschema (extra
`validate`). Tools used only to build or test (esbuild, Playwright, pytest, ruff, PyYAML, hatchling) are not
shipped either.

## Licence texts

### preact

```
The MIT License (MIT)

Copyright (c) 2015-present Jason Miller

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### @preact/signals and @preact/signals-core

```
The MIT License (MIT)

Copyright (c) 2022-present Preact Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
