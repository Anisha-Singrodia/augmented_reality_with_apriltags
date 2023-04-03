import numpy as np

def est_pixel_world(pixels, R_wc, t_wc, K):
    """
    Estimate the world coordinates of a point given a set of pixel coordinates.
    The points are assumed to lie on the x-y plane in the world.
    Input:
        pixels: N x 2 coordiantes of pixels
        R_wc: (3, 3) Rotation of camera in world
        t_wc: (3, ) translation from world to camera
        K: 3 x 3 camara intrinsics
    Returns:
        Pw: N x 3 points, the world coordinates of pixels
    """

    ##### CODE START #####
    R_wc = np.linalg.inv(R_wc)
    t_wc = -np.matmul(R_wc, t_wc).reshape(3,1)
    Rt = np.concatenate([R_wc[:,0:-1],t_wc],axis=1)
    H = np.matmul(K,Rt)
    H_inv = np.linalg.inv(H)
    pix_homo = np.concatenate([pixels[:,:],np.ones((pixels.shape[0],1))],axis=1)
    Pw = np.zeros((pixels.shape[0], 3))
    for i in range(pixels.shape[0]):
        Pw[i,:] = np.transpose(np.matmul(H_inv,np.transpose(pix_homo[i,:])))
        Pw[i,:] = Pw[i,:]/Pw[i,-1]
        Pw[i,-1] = 0
    ##### CODE END #####
    return Pw
