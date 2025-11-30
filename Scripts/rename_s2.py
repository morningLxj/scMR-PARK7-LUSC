
import os
from pathlib import Path
import shutil

BASE_DIR = Path("My_MR_Project/FigureS2")
MAPPING = {
    "FigureS2_PARK7_Sensitivity_Panel": "FigureS2A_PARK7_Sensitivity",
    "FigureS2_CTSW_Sensitivity_Panel": "FigureS2B_CTSW_Sensitivity",
    "FigureS2_TMEM50A_Sensitivity_Panel": "FigureS2C_TMEM50A_Sensitivity"
}

for old_stem, new_stem in MAPPING.items():
    for ext in ['.png', '.pdf']:
        old_file = BASE_DIR / (old_stem + ext)
        new_file = BASE_DIR / (new_stem + ext)
        if old_file.exists():
            print(f"Renaming {old_file.name} -> {new_file.name}")
            shutil.move(old_file, new_file)
        else:
            print(f"Warning: {old_file.name} not found")

# Also copy to project root or Manuscript folder if needed?
# Usually separate folder is fine.
