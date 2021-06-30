# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:35:27 2021

@author: dingxu
"""

import ftplib 

site_address = 'https://archive.stsci.edu/missions/tess/tid/s0001/0000/0000/'

ftp = ftplib.FTP(source_address = site_address)

print(ftp.getwelcome())
print(ftp.nlst())