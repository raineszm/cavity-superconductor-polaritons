from invoke import task


@task
def dropdocs(c):
    with c.cd("python/doc/"):
        c.run("make latexpdf")
        c.run("cp build/latex/*.pdf ~/Dropbox/Research/Collab/Cavity/Notes/BSNumerics/")

