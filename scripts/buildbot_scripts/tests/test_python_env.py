#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pkgutil
li = list(pkgutil.iter_modules())
print(li)

import pandas as pd

data = pd.DataFrame()

print("successfully instantiated a pandas DataFrame!")
