#!/usr/bin/env python

import click

import path_helper
import runner.runners
from runner.options import default_options


@click.group()
def main():
    pass


@main.command()
@default_options
@click.option('--parts/--no-parts', default=False)
@click.option('--gl/--no-gl', default=True)
def spinhall(fname, **kwargs):
    runner.runners.HallRunner.run(fname,
                                **kwargs)


@main.command()
@default_options
def edelstein(fname, **kwargs):
    runner.runners.EERunner.run(fname, **kwargs)


@main.command()
@default_options
def susceptibility(fname, **kwargs):
    runner.runners.SuscRunner.run(fname, **kwargs)


@main.command()
@default_options
def ginzburglandau(fname, **kwargs):
    runner.runners.GinzburgLandauRunner.run(fname, **kwargs)


if __name__ == '__main__':
    main()
