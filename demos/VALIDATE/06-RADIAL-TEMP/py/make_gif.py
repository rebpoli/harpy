#!/usr/bin/env -S python

import imageio.v2 as imageio
from PIL import Image

import glob
import os, re


dname = "sigtotxx_and_temp"
files =  glob.glob ( f"{dname}/*.png" )
natsort = lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]
files = sorted(files, key=natsort)

frames = []
for file in files:
    img = Image.open(file).convert('RGBA')  # Ensure RGBA for better color handling

    # Convert to palette mode with adaptive palette (max 256 colors), with dithering
    img = img.convert('P', palette=Image.ADAPTIVE, colors=256)

    frames.append(img)


# Open all images
# images = [Image.open(img) for img in files]

# Save as GIF
frames[0].save(
    f"{dname}.gif",
    save_all=True,
    append_images=frames[1:],
    duration=200,  # duration per frame in ms
    loop=1         # loop=0 means infinite loop
)
