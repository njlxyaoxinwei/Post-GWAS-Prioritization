import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"includes": ["scipy.special._ufuncs_cxx"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
# if sys.platform == "win32":
#     base = "Win32GUI"

setup(  name = "algo",
        version = "0.1",
        description = "Algorithm",
        options = {"build_exe": build_exe_options},
        executables = [Executable("algo.py", base=base)])
