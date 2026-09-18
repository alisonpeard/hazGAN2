!#/bin/zsh

git filter-repo --force --blob-callback '
import json
try:
    nb = json.loads(blob.data)
except Exception:
    nb = None
if nb and "cells" in nb:
    changed = False
    for cell in nb["cells"]:
        if cell.get("cell_type") == "code":
            if cell.get("outputs"):
                cell["outputs"] = []
                changed = True
            if cell.get("execution_count") is not None:
                cell["execution_count"] = None
                changed = True
    if changed:
        blob.data = json.dumps(nb, indent=1).encode()
'
