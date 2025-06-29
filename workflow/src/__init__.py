"""
>>> from src import funcs;print(funcs.identity)
>>> export SNAKEMAKE_PROJECT=bayofbengal_era5
>>> python -c "from src import funcs;print(funcs.identity);print(funcs.wind_speed);"
"""
import os
import sys
from pathlib import Path
from .python import funcs
from .python import datasets


def load_project_funcs(project_name=None):
    """Load project-specific functions into the funcs module"""

    if project_name is None:
        project_name = os.environ.get("SNAKEMAKE_PROJECT")
    
    if not project_name:
        return
    
    project_path = Path(__file__).parents[2] / "projects" / project_name / "src"

    if project_path.exists():
        sys.path.insert(0, str(project_path))
        try:
            import funcs as project_funcs
            
            # Add project functions to the funcs module we imported above
            for name in dir(project_funcs):
                if not name.startswith('_'):
                    setattr(funcs, name, getattr(project_funcs, name))
        except ImportError:
            pass
        finally:
            sys.path.remove(str(project_path))


load_project_funcs()