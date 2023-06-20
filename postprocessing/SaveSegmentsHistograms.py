# %% Import zeoPy tool
from zeoPy import zeoPy

# %% Import required packages

import pandas as pd
import numpy as np
import os 

# %%
f = zeoPy()
#uniqueSegmentsIZA,bounds = f.getUniqueSegments("../MaterialsInfo",removeIsolated=True)
f.computeMaterialsHistograms("../MaterialsInfo", "../IZASegmentHistograms/",150,removeIsolated=True)