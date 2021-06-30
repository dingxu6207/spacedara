# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astroquery.mast import Catalogs


target_name = "2MASS J04215943+1932063"
search_radius_deg = 0.2

# Query the TESS Input Catalog centered on HD 209458 with a 0.2 degree radius.
catalogTIC = Catalogs.query_object(target_name, radius=search_radius_deg, catalog="TIC")

# Print out the number of returned rows.
print("Number of TIC objects within %f deg of %s: %u" % (search_radius_deg, target_name, len(catalogTIC)))

print(type(catalogTIC))