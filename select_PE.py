"""
Print list of PE numbers whose fort.14 files contain the point x, y
in lat/lon.

Example:
    `select_PE.py 120m_mesh/ x y`
"""

#!/usr/bin/env python3
import utils_beach as ub
import os
import pandas as pd

if __name__ == "main":
    y = 8.782
    x = -81.304

    proj = "120m_nolevee"
    folder = os.path.abspath(proj)
    domains = list(filter(lambda x: "PE" in x, os.listdir(folder)))


    for domain in domains:
        f14 = ub.fort14(os.path.join(folder, domain + "/fort.14"))
        if f14.is_in_mesh(x, y):
            print(domain)
