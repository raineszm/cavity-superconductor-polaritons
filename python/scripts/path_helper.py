import pathlib
import sys

lib_dir = pathlib.Path(__file__).absolute().parents[1] / 'lib'

sys.path.append(str(lib_dir))
