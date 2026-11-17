# Third-party notices

The scratch directory ownership and bounded external merge-sort design in
`src/store.hpp`, and parts of the build/bootstrap conventions, are adapted
from [gwaslocusfinder](https://github.com/bgorm/gwaslocusfinder).

Copyright (c) 2026 gwaslocusfinder contributors

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

HTSlib and bcftools are external MIT/Expat-licensed dependencies from
[samtools](https://www.htslib.org/). The plugin uses their headers and resolves
HTSlib symbols from its host. It does not bundle or link an additional HTSlib.
The optional bootstrap downloads the HTSlib source archive, which retains its
own license notices.

The small gene annotation fixture under `tests/fixtures/` contains unmodified
Ensembl release 116 records. Source URLs, checksums, selection rules, and
coverage are documented in [the fixture notes](tests/fixtures/README.md).
