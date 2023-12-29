import matplotlib as mpl
import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, MultiPoint, LineString, MultiLineString
import shapely.plotting
import shapely.affinity as affinity
import shapely

import geopandas as gpd

import numpy as np

import math
import collections

from .geom_helpers import *

class BandSlicer: #TODO: figure out inheritance model for parent slicer class
    def __init__(self, geom, tolerance, n_bands, band_width, piece_width):
        self.tolerance = tolerance #TODO: standardize this somehow
        self.n_bands = n_bands       #TODO: decide this based on target piece num.
        self.band_width = band_width #TODO: decide this based on bands
        self.piece_width = piece_width       #TODO: decide this based on target piece num.

        self.geom = self.preprocess(geom)

        self.rng = np.random.default_rng()
        self.longitudes = self._longitudes()
        self.band_cuts  = self.cut_bands()

    def _longitudes(self):
        buffer_vals = np.arange(-self.band_width, -self.band_width*self.n_bands, -self.band_width)
        return [self.geom.boundary] + [self.geom.buffer(val).boundary for val in buffer_vals]

    def preprocess(self, geom: shapely.Polygon):
        if isinstance(geom, gpd.GeoDataFrame):
          geom = geom.loc[0, "geometry"]
        elif isinstance(geom, gpd.GeoSeries):
          geom = geom[0]
        
        geom = self.normalize_geom(geom)
        geom = geom.simplify(self.tolerance)
        #geom = geom.buffer(-self.band_width/5)

        return geom

    def plot(self):
        _, ax = plt.subplots()
        for long in self.longitudes:
            #long.plot(ax=ax, facecolor="none")
            shapely.plotting.plot_line(long, ax=ax,  linewidth=1, add_points=False)

        for band in self.band_cuts:
            shapely.plotting.plot_line(MultiLineString(band), linewidth=1, ax=ax, add_points = False)

        ax.axis("equal")
        ax.set_axis_off()

    def cut_band(self, inner_lng_idx: int, n_cuts: int, offset: int):
        cuts = []
        used_pts = []

        inner_long = self.longitudes[inner_lng_idx]
        outer_long = self.longitudes[inner_lng_idx - 1]

        ideal_separation = self.piece_width

        for i, coord in enumerate(inner_long.coords):
            edges = oriented_edges(inner_long, i)
            int_angle = interior_angle(inner_long, i)

            if int_angle < 120:
                cuts += [perp_line(coord, edges[0], self.band_width), perp_line(coord, edges[1], self.band_width)]
                used_pts += [coord]

            elif int_angle < 160:
                test_line = perp_line(coord, edges[0]/norm(edges[0]) + edges[1]/norm(edges[1]), 10*self.band_width)
                new_end_pt = shapely.intersection(test_line, outer_long).representative_point()
                
                cuts += [LineString([coord, new_end_pt])]
                used_pts += [coord]
        
        remaining_cuts = n_cuts - len(cuts)

        if remaining_cuts == n_cuts:
            distances = np.linspace(0, inner_long.length, num=n_cuts, endpoint=False)
            points = [inner_long.interpolate(distance) for distance in distances]
            cuts += [shapely.shortest_line(pt, outer_long) for pt in points]

            return cuts

        distance_breakpoints = [inner_long.line_locate_point(Point(pt)) for pt in used_pts]

        for cur, next in zip(distance_breakpoints[:-1], distance_breakpoints[1:-1] + [inner_long.length]):
            if abs(next - cur) < 1.5*ideal_separation:
              continue
                   
            n = math.ceil(abs((next - cur))/ideal_separation) + 1
            segment_distances = np.linspace(cur, next, num=n)[1:-1]
            points = [inner_long.interpolate(distance) for distance in segment_distances]
            cuts += [shapely.shortest_line(pt, outer_long) for pt in points]
        
        print(n_cuts - len(cuts))

        return cuts

    def cut_bands(self):
        outside_spacing = self.piece_width
        band_idxs = range(0,self.n_bands-1)

        band_cut_res = [
            (idx+1, math.ceil(self.longitudes[idx+1].length/outside_spacing))
            for idx in band_idxs]

        return [self.cut_band(*band_data, 0) for band_data in band_cut_res]

    @staticmethod
    def normalize_geom(geom: shapely.Polygon):
        xmin, ymin, xmax, ymax = geom.bounds
        max_dim = max(xmax-xmin, ymax-ymin)

        transed = affinity.translate(geom, xoff = -xmin, yoff = -ymin)
        return affinity.scale(transed, xfact = 1/max_dim, yfact = 1/max_dim, origin=(0,0))