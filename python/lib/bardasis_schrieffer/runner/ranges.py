import click
import numpy as np


class Range(click.ParamType):

    name = "Range"

    def convert(self, value, param, ctx):
        if isinstance(value, np.ndarray):
            return value
        split = value.split(":")
        if len(split) != 3:
            self.fail("Incorrectly formatted range {}".format(value), param, ctx)
        lo, hi, n = split
        return np.linspace(float(lo), float(hi), int(n))
