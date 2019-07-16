from pathlib import Path

from doit.action import CmdAction

NOTES_DIR = Path("~/Dropbox/Research/Collab/Cavity/Notes/").expanduser()


def task_dropdocs():
    return {
        "actions": [
            CmdAction(cmd, cwd="python/doc")
            for cmd in [
                "make latexpdf html",
                f'cp build/latex/*.pdf {NOTES_DIR / "BSNumerics"}',
                f'cp -R build/html {NOTES_DIR / "BSNumerics"}',
            ]
        ]
    }


def task_pull_notes():
    texfile = "Bardasis-Schrieffer Polaritons.pdf"
    return {
        "actions": [
            CmdAction(
                f'cp "{NOTES_DIR / "Bardasis-Schrieffer Polaritons" / texfile}" .',
                cwd="python/doc",
            )
        ]
    }
