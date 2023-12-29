import matplotlib.pyplot as plt

from shapely.geometry import Point, Polygon, MultiPoint, LineString, MultiLineString
import shapely.plotting
import shapely.affinity as affinity
import shapely

import geopandas as gpd

import numpy as np

import math

import warnings

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

    def preprocess(self, geom: shapely.Polygon):
        if isinstance(geom, gpd.GeoDataFrame):
          geom = geom.loc[0, "geometry"]
        elif isinstance(geom, gpd.GeoSeries):
          geom = geom[0]
        
        geom = self.normalize_geom(geom)
        geom = geom.simplify(self.tolerance)

        return geom

    def _longitudes(self):
        buffer_vals = np.arange(-self.band_width, -self.band_width*self.n_bands, -self.band_width)
        long_list = [self.geom.boundary] + [self.geom.buffer(val).boundary for val in buffer_vals]

        if any([geom.is_empty for geom in long_list]):
            long_list = [x for x in long_list if not x.is_empty]
            warnings.warn("Too many longitudes, empty longitudes discarded")
        
        return long_list

    def cut_bands(self):
        outside_spacing = self.piece_width
        band_idxs = range(0,self.n_bands-1)

        band_cut_res = [
            (idx+1, math.ceil(self.longitudes[idx+1].length/outside_spacing))
            for idx in band_idxs] ## List of tuples of golw (longitude index, # of cuts that fit to give desired piece width)

        return [self.cut_band(*band_data, 0) for band_data in band_cut_res]


    def cut_band(self, inner_lng_idx: int, n_cuts: int, offset=0): ## NOTE: offset not currently being used
        cuts = []
        used_pts = []

        inner_long = self.longitudes[inner_lng_idx]
        outer_long = self.longitudes[inner_lng_idx - 1]

        if isinstance(outer_long, shapely.MultiLineString):
            nearby_geoms = list(outer_long.geoms)
        else:
            nearby_geoms = [outer_long]

        if isinstance(inner_long, shapely.MultiLineString):
            print("oopises, i should implement this :)")
        else:
            for cut in self.identify_cut_points(inner_long, self.piece_width):
                cuts += self.line_at_cut(cut, inner_long, nearby_geoms)
            
        return cuts
        
    def plot(self):
        _, ax = plt.subplots()
        for long in self.longitudes:
            shapely.plotting.plot_line(long, ax=ax,  linewidth=1, add_points=False)

        for band in self.band_cuts:
            shapely.plotting.plot_line(MultiLineString(band), linewidth=1, ax=ax, add_points = False)

        ax.axis("equal")
        ax.set_axis_off()

    @staticmethod
    def normalize_geom(geom: shapely.Polygon): ## TODO: This should either move to geom_helper, or the geom_helper methods should move here
        xmin, ymin, xmax, ymax = geom.bounds
        max_dim = max(xmax-xmin, ymax-ymin)

        transed = affinity.translate(geom, xoff = -xmin, yoff = -ymin)
        return affinity.scale(transed, xfact = 1/max_dim, yfact = 1/max_dim, origin=(0,0))
    
    @staticmethod
    def identify_cut_points(geom: shapely.LineString, separation: float):
        cut_points = []

        for i, coord in enumerate(geom.coords): ## start by cutting at vertices
            pt = {"distance": geom.line_locate_point(Point(coord)), "coords": coord, "type": "direction"}

            edges = oriented_edges(geom, i)
            int_angle = interior_angle(geom, i)

            if int_angle < 120:
                pt["directions"] = [perp_vec(edges[0]), perp_vec(edges[1])]

            elif int_angle < 160: 
                pt["directions"] = [perp_vec(edges[0]/norm(edges[0]) + edges[1]/norm(edges[1]))]
            
            else: continue

            cut_points += [pt]
        
        if len(cut_points) != 0:
            distance_breakpoints = [pt["distance"] for pt in cut_points]
        else:
            distance_breakpoints = [0, geom.length]

        for cur, next in zip(distance_breakpoints[:-1], distance_breakpoints[1:-1] + [geom.length]): ## then cut up the segments between verts
            if abs(next - cur) < 2*separation:
              continue
                   
            n = math.ceil(abs((next - cur))/separation) + 1
            segment_distances = np.linspace(cur, next, num=n)[1:]
            cut_points += [{
                "distance": d,
                "coords":   geom.interpolate(d).coords,
                "type": "nearest"
            } for d in segment_distances]
        
        return cut_points

    @classmethod
    def line_at_cut(cls, cut, geom: shapely.LineString, nearby_geoms: list[shapely.LineString] = []):
        if cut["type"] == "direction":
            test_lines = [dir_line(cut["coords"], dir, math.sqrt(2)) for dir in cut["directions"]]
            return [cls.clip_line(line, nearby_geoms) for line in test_lines]
        elif cut["type"] == "nearest":
            test_lines = [shapely.shortest_line(Point(cut["coords"]), geom) for geom in nearby_geoms]
            return [min(test_lines, key = lambda line: line.length)]

    @staticmethod
    def clip_line(line: shapely.LineString, nearby_geoms: list[shapely.LineString]):
        base_pt = Point(line.coords[0])

        intersections = [shapely.intersection(line, geom) for geom in nearby_geoms]

        if not intersections:
            return line
        
        closest = min(intersections, key=lambda x: shapely.shortest_line(base_pt, x).length)
        return shapely.shortest_line(base_pt, closest)
