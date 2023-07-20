import glob
import os
import shutil

from labw_utils.devutils import myst_nb_helper

os.environ['SPHINX_BUILD'] = '1'  # Disable chronolog and others.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)


def scan_dir(path_to_scan: str):
    """
    Recursively list files in directories.

    <https://www.sethserver.com/python/recursively-list-files.html>
    """

    files = []
    dirlist = [path_to_scan]
    while len(dirlist) > 0:
        for (dirpath, dirnames, filenames) in os.walk(dirlist.pop()):
            dirlist.extend(dirnames)
            files.extend(map(lambda n: os.path.join(*n), zip([dirpath] * len(filenames), filenames)))
    return files


def copy_doc_files(from_path: str, to_path: str):
    """
    Copy items to project root
    """
    os.makedirs(to_path, exist_ok=True)
    for name in glob.glob(from_path):
        shutil.copy(name, to_path + os.sep)


copy_doc_files(os.path.join(ROOT_DIR, '*.md'), os.path.join(THIS_DIR, "_root"))

myst_nb_helper.convert_ipynb_to_myst(
    THIS_DIR,
    hooks=[myst_nb_helper.shell_filter]
)

myst_nb_helper.generate_cli_docs(
    os.path.join(THIS_DIR, "cli_docs.toml"),
    os.path.join(THIS_DIR, "_cli_docs")
)
