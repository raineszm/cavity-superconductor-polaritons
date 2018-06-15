from invoke import task
from pathlib import Path

NOTES_DIR = Path("~/Dropbox/Research/Collab/Cavity/Notes/").expanduser()


@task
def dropdocs(c):
    with c.cd("python/doc/"):
        c.run("make latexpdf html")
        c.run(f'cp build/latex/*.pdf {NOTES_DIR / "BSNumerics"}')
        c.run(f'cp -R build/html {NOTES_DIR / "BSNumerics"}')


@task
def pull_notes(c):
    with c.cd("python/doc/source/ext/"):
        texfile = "Bardasis-Schrieffer Polaritons.pdf"
        c.run(f'cp "{NOTES_DIR / "Bardasis-Schrieffer Polaritons" / texfile}" .')
