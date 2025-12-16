import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from part2d.part2d import part2d

def clear_files():
    for f in glob.glob("part2b/nu0.01/scheme1/*"):
        os.remove(f)
    for f in glob.glob("part2b/nu0.01/scheme2/*"):
        os.remove(f)
    for f in glob.glob("part2b/nu1.0/scheme1/*"):
        os.remove(f)
    for f in glob.glob("part2b/nu1.0/scheme2/*"):
        os.remove(f)
    

clear_files()
part2d()


