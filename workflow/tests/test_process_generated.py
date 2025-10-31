# %%
"""WIP: tests for process_generated.py"""
import sys
from pathlib import Path

workflow_dir = Path(__file__).parent.parent
scripts_dir = workflow_dir / "scripts"

sys.path.insert(0, str(workflow_dir))
sys.path.insert(0, str(scripts_dir))

import process_generated
# %%
import numpy as np
import xarray as xr

mri = 1000
p = 1 - (1 / mri)
# 
fields = {}.
ref = xr.open_dataset(input.training_data)
x_ref = ref["anomaly"].values
# get dimensions of ref
h, w, c = x_ref.shape[1], x_ref.shape[2], x_ref.shape[3]
u = p * np.ones((1, h, w, c))  # uniform array
process_generated.transform_to_physical(u, ref, False)
# %%
