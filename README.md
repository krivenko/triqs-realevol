realevol - Real time evolution solver based on TRIQS
====================================================

Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek

Installation
------------

- Install the [TRIQS](http://triqs.github.io/triqs/3.0.x/install.html) library version 3.0.x.
- `source <path_to_triqs_install_dir>/share/triqsvars.sh`
- `git clone https://github.com/krivenko/triqs-realevol.git triqs-realevol.git`
- `mkdir triqs-realevol.build && cd triqs-realevol.build`
- `cmake ../triqs-realevol.git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<path_to_install_dir>`
- `make`
- `make test`
- `make install`
