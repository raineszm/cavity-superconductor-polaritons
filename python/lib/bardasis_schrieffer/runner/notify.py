import pushbullet
import click
import pendulum
import textwrap
import pathlib


def _get_pb():
    key_file = pathlib.Path("~/.pushbullet.key").expanduser()
    if not key_file.exists():
        click.secho("Unable to find ~/.pushbullet.key", fg="red", err=True)
        return

    return pushbullet.Pushbullet(key_file.read_text().strip())


def notify_success(fname, N):
    pb = _get_pb()

    if not pb:
        return

    pb.push_note(
        "Bardasis Schrieffer run completed",
        textwrap.dedent(
            f"""\
                 {fname} completed {pendulum.now().format('LLLL')}

                 Generated {N} datapoints
                 """
        ),
    )


def notify_failure(fname, Nfail, exc):
    pb = _get_pb()

    if not pb:
        return

    pb.push_note(
        "Bardasis Schrieffer run failed",
        textwrap.dedent(
            f"""\
                 {fname} failed {pendulum.now().format('LLLL')}

                 with {Nfail} errors

                 Final error: {exc}
                 """
        ),
    )

