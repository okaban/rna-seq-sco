"""
Thin wrapper around 00_shared_utils.py for use in new figure scripts.
Imports shared utils functions/constants without duplicating them.
"""
from pathlib import Path
import sys

# Add this directory to path so 00_shared_utils is importable
_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

# Import everything from 00_shared_utils
from importlib import import_module
_mod = import_module('00_shared_utils')

apply_style = _mod.apply_style
add_panel_label = _mod.add_panel_label
save_figure = _mod.save_figure
FIG_SUP_DIR = _mod.FIG_SUP_DIR
FIG_DIR = _mod.FIG_DIR

# Re-export color constants
COL_4mC = _mod.COL_4mC
COL_6mA = _mod.COL_6mA
COL_BOTH = _mod.COL_BOTH
COL_GRAY = _mod.COL_GRAY
COL_DARK = _mod.COL_DARK
COL_SHIELDED = _mod.COL_SHIELDED
COL_EXPOSED = _mod.COL_EXPOSED
COL_ACTIVATION = _mod.COL_ACTIVATION
COL_REPRESSION = _mod.COL_REPRESSION
TP_LABELS = _mod.TP_LABELS
