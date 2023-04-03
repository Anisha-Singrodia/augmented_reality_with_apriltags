from est_homography import est_homography
import numpy as np

def PnP(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-N-Point problem with collineation assumption, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3, ) numpy array describing camera translation in the world (t_wc)

    """

    ##### CODE START #####

    # Homography Approach
    H = est_homography(Pw[:,0:2], Pc)
    z = H[2][2]
    H = H/z
    K_inv_H = np.matmul(np.linalg.inv(K), H)
    # Following slides: Pose from Projective Transformation
    H_dash = np.zeros((3,3))
    H_dash[:, 0] = K_inv_H[:, 0]
    H_dash[:, 1] = K_inv_H[:, 1]
    H_dash[:, 2] = np.cross(np.transpose(K_inv_H[:, 0]), np.transpose(K_inv_H[:,1]))
    [U, S , Vt ] = np.linalg.svd(H_dash)
    mat = np.eye(3)
    mat[2][2] = np.linalg.det(np.matmul(U, Vt))
    R = np.transpose(np.matmul(np.matmul(U, mat), Vt))
    t = -np.matmul(R, (K_inv_H[:,2]/(np.linalg.norm(K_inv_H[:,0]))))

    ##### CODE END #####
    return R, t               
