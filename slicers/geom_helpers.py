import shapely
import numpy as np
from numpy.linalg import norm
import math

def oriented_edges(poly: shapely.Geometry, pt_idx: int):
  if shapely.get_dimensions(poly) == 1:
      coords = poly.coords
  elif shapely.get_dimensions(poly) == 2:
      coords = poly.boundary.coords
  else:
      raise ValueError(f"This method does not support {shapely.get_dimensions(poly)}d geometry types")

  pt_idx = pt_idx % (len(coords) - 1)

  v2 = np.array(coords[pt_idx + 1]) - np.array(coords[pt_idx])

  if pt_idx == 0:
    v1 = np.array(coords[pt_idx]) - np.array(coords[pt_idx - 2])
  else:
    v1 = np.array(coords[pt_idx]) - np.array(coords[pt_idx - 1])

  return (v1, v2)

def interior_angle(poly: shapely.Geometry, pt_idx: int):
  edges = oriented_edges(poly, pt_idx)
  angle = math.degrees(-np.arctan2(-edges[0][1], -edges[0][0]) + np.arctan2(edges[1][1],edges[1][0]))

  if angle < 0:
    angle += 360
  
  return angle

def perp_vec(vec: list[float]):
  return np.array([-vec[1], vec[0]])

def dir_line(base_pt: list[float], dir_vec: list[float], length: float):
  base_pt, vector = (np.array(base_pt), np.array(dir_vec))
  unit     = dir_vec/norm(dir_vec)

  end_pt = base_pt + length*unit
  return shapely.LineString([base_pt, end_pt])

def perp_line(base_pt: list[float], vector: list[float], length: float):
  base_pt, vector = (np.array(base_pt), np.array(vector))

  unit     = vector/norm(vector)
  perp_vec = np.array([-unit[1], unit[0]])

  end_pt = base_pt + length*perp_vec
  return shapely.LineString([base_pt, end_pt])

#TODO: force pt and vec inputs to be tuples and or np arrays