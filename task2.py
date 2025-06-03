import importlib
import sys, os
import pathlib

root = pathlib.Path("C:\\Program Files\\Lumerical")

# Recursive search for a python/lumapi.py file in the specified directory
for dirpath, dirnames, filenames in os.walk(root):
    if 'lumapi.py' in filenames:
        sys.path.append(dirpath)
        print("Found lumapi at:", os.path.join(dirpath, 'lumapi.py'))
        break
else:
    raise ImportError("lumapi.py not found in the specified directory tree.")

import lumapi
print(lumapi.__file__)