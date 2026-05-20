# realevol - Real time evolution solver based on TRIQS

Copyright (C) 2014-2026, I. Krivenko, M. Danilov, P. Kubiczek

## Installation

### `conda`/`mamba` package from Anaconda.org

```
conda install -c conda-forge krivenko::triqs_realevol
mamba install -c conda-forge krivenko::triqs_realevol
```

### Compilation from sources

- Install the [TRIQS](http://triqs.github.io/triqs/3.3.x/install.html) library version 3.3.x.
- Make sure you have [ARPACK-NG](https://github.com/opencollab/arpack-ng) version 3.8.0 or newer installed somewhere.
- `source <path_to_triqs_install_dir>/share/triqsvars.sh`
- `git clone https://github.com/krivenko/triqs-realevol.git triqs-realevol.git`
- `mkdir triqs-realevol.build && cd triqs-realevol.build`
- `cmake ../triqs-realevol.git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<path_to_install_dir> -Darpack-ng_DIR=<ARPACK-NG_installation_prefix>/lib/cmake`
- `make`
- `make test`
- `make install`

## Acknowledgements

Development of this software was partially funded by the European Research
Council (ERC) under the European Union's Horizon 2020 research and innovation
programme (Grant Agreement No.
[854843-FASTCORR](https://cordis.europa.eu/project/id/854843)).

## License

realevol is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

realevol is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
realevol (in the file LICENSE.txt in this directory). If not, see
<http://www.gnu.org/licenses/>.
