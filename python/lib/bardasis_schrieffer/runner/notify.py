import pushbullet
import click
import pendulum
import textwrap
import pathlib


def push_notification(fname, N):
    key_file = pathlib.Path("~/.pushbullet.key").expanduser()
    if not key_file.exists():
        click.secho("Unable to find ~/.pushbullet.key", fg="red", err=True)
        return

    pb = pushbullet.Pushbullet(key_file.read_text().strip())

    pb.push_note(
        "Bardasis Schrieffer run completed",
        textwrap.dedent(
            f"""\
                 {fname} completed {pendulum.now().format('LLLL')}

                 Generated {N} datapoints
                 """
        ),
    )

