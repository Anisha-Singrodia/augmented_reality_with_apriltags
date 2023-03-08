import numpy as np

def est_homography(X, Y):
    """
    Calculates the homography H of two planes such that Y ~ H*X
    If you want to use this function for hw5, you need to figure out
    what X and Y should be.
    Input:
        X: 4x2 matrix of (x,y) coordinates
        Y: 4x2 matrix of (x,y) coordinates
    Returns:
        H: 3x3 homogeneours transformation matrix s.t. Y ~ H*X

    """

    ##### STUDENT CODE START #####
    A = np.zeros([8, 9])
    for i in range(4):
        A[2*i][:] = [-X[i][0], -X[i][1], -1, 0, 0, 0, X[i][0]*Y[i][0], X[i][1]*Y[i][0], Y[i][0]]
        A[2*i+1][:] = [0, 0, 0, -X[i][0], -X[i][1], -1, X[i][0]*Y[i][1], X[i][1]*Y[i][1], Y[i][1]]
    [U, S , Vt ] = np.linalg.svd(A)
    H = Vt[:][-1]
    H = H.reshape(3,3)

    ##### STUDENT CODE END #####

    return H
