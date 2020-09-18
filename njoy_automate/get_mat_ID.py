#!/usr/bin/env python

import sys
import numpy as np

file = sys.argv[1]
with open(file, "r") as f:
    for line in f:
        if "MATERIAL" in line:
            tokens = np.array(line.strip().split())
            ind = np.argwhere(tokens == "MATERIAL").item() + 1
            mat = tokens[ind]
            print(mat)
            break
sys.exit(0)
