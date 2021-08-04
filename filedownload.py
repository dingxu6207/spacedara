# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 17:36:56 2021

@author: dingxu
"""

import os
path = 'I:\\TESSDATA\\section1\\'
dirpath = 'I:\\TESSDATA\\variable\\section1\\'
for root, dirs, files in os.walk(path):
   for file in files:
       strfile = os.path.join(root, file)
       if (strfile[-3:] == '.sh'):
           os.system(strfile)