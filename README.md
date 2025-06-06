# SALMON software document

SALMON software document powered by Sphinx documentation builder.

SALMON: http://salmon-tddft.jp/

Sphinx: http://www.sphinx-doc.org/en/master/

## Setup

Please install the Sphinx and functional packages with `pip` package manager.

    $ pip install sphinx sphinx_rtd_theme sphinxcontrib-bibtex

We recommend the `sphinx-autobuild` package.

    $ pip install sphinx-autobuild
    $ sphinx-autobuild source build/html
    $ open http://127.0.0.1:8000/

## Build

In Linux and Mac OS X, you can use the GNU make for building the document.
If you use Windows, `make.bat` is assisted to build the document.

## License

SALMON is available under Apache License version 2.0.

    Copyright 2017-2025 SALMON developers
    
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
  
       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
