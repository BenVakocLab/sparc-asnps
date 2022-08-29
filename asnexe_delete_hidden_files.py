# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:22:59 2022

@author: vlab_eye
"""

import sys, os
import numpy as np
from matplotlib import pyplot as plt

path = 'E:\\ps-artifact'

list_files = []

for r, d, f in os.walk(path):
    for file in f:
        if '.DS_Store' in file:
            path_to_del = os.path.join(r, file)
            list_files.append(path_to_del)
            os.remove(path_to_del)
            
if len(list_files)<1:
    print("nothing found")
else:
    for f in list_files:
        print(f)