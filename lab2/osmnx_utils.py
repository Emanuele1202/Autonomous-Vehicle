import numpy as np
import numpy as np


def G2C_RGF(longitude, latitude, h) :

  GRS_a = 6378137
  GRS_f = 1/298.257222101
  GRS_b = GRS_a*(1-GRS_f)
  GRS_e = np.sqrt((np.power(GRS_a,2) - np.power(GRS_b,2)) / np.power(GRS_a,2))
  N = GRS_a / np.sqrt(1 - np.power(GRS_e,2) * np.power(np.sin(latitude),2))

  # TODO : compute the coordinates of the point in the RGF93 frame where RGF93 is a local frame
  
  X = (N + h) * np.cos(latitude) * np.cos(longitude)
  Y = (N + h) * np.cos(latitude) * np.sin(longitude)
  Z = (N * (1 - np.power(GRS_e,2)) + h) * np.sin(latitude)


  return np.array([X, Y, Z])

def C_RGF2ENU(P_ECEF, l, phi, h) :
  eO = G2C_RGF(l, phi, h)
  P_ECEF = P_ECEF[np.newaxis].T # (3,1)

  # TODO : compute the coordinates of the point in the local frame
  # i do it with by hand , each term comes from the matrix multiplication between the rotation matrix and the vector
  oAe = np.eye(3)
  oAe[0,0] = -np.sin(l)
  oAe[0,1] = np.cos(l)
  oAe[1,0] = -np.sin(phi)*np.cos(l)
  oAe[1,1] = -np.sin(phi)*np.sin(l)
  oAe[1,2] = np.cos(phi)
  oAe[2,0] = np.cos(phi)*np.cos(l)
  oAe[2,1] = np.cos(phi)*np.sin(l)
  oAe[2,2] = np.sin(phi)

  P_ENU = oAe @ (P_ECEF -  eO[np.newaxis].T)

  return P_ENU[:, 0]



def gdfs_to_local(nodes_proj, path):
  ''' This function takes a list of nodes and returns the coordinates of the nodes in a local frame
  Args:
    nodes_proj (gpd.GeoDataFrame): a GeoDataFrame of nodes
    path (list): a list of nodes' ids
  Returns:
    waypoints_ENU (np.array): an array of shape (len(path), 3) containing the coordinates of the nodes in the local frame
  '''
  l = nodes_proj.at[nodes_proj.index[0],"x"] * np.pi / 180
  phi = nodes_proj.at[nodes_proj.index[0],"y"] * np.pi / 180
  h0 = 342.

  # put the points in a frame in meter :
  waypoints_ECEF = np.zeros((nodes_proj.loc[path]["x"].size, 3))
  waypoints_ENU = np.zeros((nodes_proj.loc[path]["x"].size, 3))
  
  for i in range(nodes_proj.loc[path]["x"].size):
      longitude = nodes_proj.loc[path[i]]["x"] * np.pi / 180
      latitude =  nodes_proj.loc[path[i]]["y"] * np.pi / 180

      # TODO : compute the coordinates of the point in the RGF frame
      waypoints_ECEF[i] = G2C_RGF(longitude, latitude, h0)
      # TODO : compute the coordinates of the point in the local frame
      waypoints_ENU[i] = C_RGF2ENU(waypoints_ECEF[i], l, phi, h0)

  return waypoints_ENU

