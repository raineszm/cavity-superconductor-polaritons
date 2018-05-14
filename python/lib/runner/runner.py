import multiprocessing as mp
import csv
import pathlib
import textwrap
from datetime import datetime

import numpy as np
import pushbullet
import click

import bardasis_schrieffer
from .defaults import *


def push_notification(fname, N, METADATA):
    key_file = pathlib.Path('~/.pushbullet.key').expanduser()
    if not key_file.exists():
        click.secho('Unable to find ~/.pushbullet.key', fg='red', err=True)
        return

    pb = pushbullet.Pushbullet(key_file.read_text().strip())

    pb.push_note('Bardasis Schrieffer run completed',
                 textwrap.dedent('''\
                 {fname} completed at {dtime}

                 Generate {N} datapoints
                 '''.format(fname=fname,
                            dtime=datetime.now(),
                            N=N,
                            meta=METADATA)))


def run(fname, confirm=False, **kwargs):

    if confirm:
        get_confirmation(**kwargs)

    N = 1
    ts = [0.]

    with mp.Pool() as p:
        with open(fname, 'w') as datafile:
            writer = csv.DictWriter(
                datafile,
                fieldnames=[]
            )

            writer.writeheader()

            with click.progressbar(
                    length=N,
                    fill_char=click.style('#', fg='green'),
                    empty_char=click.style('-', fg='red')) as bar:

                generator = cls.make()

                for result in p.imap_unordered(generator, ts):
                    writer.writerow(result)
                    datafile.flush()
                    bar.update(1)

    if kwargs['notify']:
        push_notification(fname, N, dict())


def show_param(name, value):
    if hasattr(value, '__len__'):
        return '{}: [{}, {}] x {}'.format(name,
                                          min(value), max(value), len(value))
    else:
        return '{}: {}'.format(name, value)


def get_confirmation(**params):
    click.secho('Please confirm parameters:', fg='green')
    for (k, v) in params.items():
        click.echo(show_param(k, v))
    click.confirm('Is this correct?', abort=True)
