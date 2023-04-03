import numpy as np

def est_Pw(s):
    """
    Estimate the world coordinates of the April tag corners, assuming the world origin
    is at the center of the tag, and that the xy plane is in the plane of the April
    tag with the z axis in the tag's facing direction. See world_setup.jpg for details.
    Input:
        s: side length of the April tag

    Returns:
        Pw: 4x3 numpy array describing the world coordinates of the April tag corners
            in the order of a, b, c, d for row order. See world_setup.jpg for details.

    """

    ##### CODE START #####
    points = np.zeros((4,3))
    points[0][:] = np.array([0,0,0])
    points[1][:] = np.array([s,0,0])
    points[2][:] = np.array([s,s,0])
    points[3][:] = np.array([0,s,0])
    print(points)
    ##### CODE END #####

    Pw = np.zeros((4,3))
    for i in range(points.shape[0]):
        for j in range(points.shape[1] - 1):
            Pw[i][j] = points[i][j] - (s/2)
    return Pw
# print(est_Pw(0.14))
