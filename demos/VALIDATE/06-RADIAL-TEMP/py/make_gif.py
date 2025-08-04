#!/usr/bin/env -S python

from PIL import Image

import glob
import os, re


files =  glob.glob ( f"sigyy_pngs/*.png" )

natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]
files = sorted(files, key=natsort)

# Open all images
images = [Image.open(img) for img in files]

# Save as GIF
images[0].save(
    "sigyy.gif",
    save_all=True,
    append_images=images[1:],
    duration=200,  # duration per frame in ms
    loop=1         # loop=0 means infinite loop
)
