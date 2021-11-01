#!/usr/bin/python
import os
import argparse
import numpy as np
import folium


class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    """Hack to allow showing default values in help and preserving line break in epilog"""

    pass


parser = argparse.ArgumentParser(
    description="Currents validation.", formatter_class=CustomFormatter
)

parser.add_argument(
    "-x", "--lons", type=float, nargs="+", help="List of longitudes. Ex: -44 -45.3 -41"
)

parser.add_argument(
    "-y", "--lats", type=float, nargs="+", help="List of latitudes: Ex: -23 -34 -28"
)

args = parser.parse_args()
print("\n")

m = folium.Map(location=[np.mean(args.lats), np.mean(args.lons)], zoom_start=2)

for lon, lat in zip(args.lons, args.lats):
    folium.Marker(location=[lat, lon]).add_to(m)

m.save("/tmp/map.html")

os.system("google-chrome /tmp/map.html")
