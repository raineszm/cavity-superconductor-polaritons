#!/usr/bin/env python

import pathlib
import sys

lib_dir = pathlib.Path(__file__).absolute().parents[1] / "lib"

sys.path.append(str(lib_dir))

import bardasis_schrieffer.runner.cli

if __name__ == '__main__':
    bardasis_schrieffer.runner.cli.main()
