# (c) 2020 by Authors
# This file is a part of the SD program.
# see LICENSE file

import subprocess


def get_git_revision_hash():
    return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()


def get_git_revision_short_hash():
    return (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )
