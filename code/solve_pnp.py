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

    ##### STUDENT CODE START #####

    # Homography Approach
    # print(Pc)
    # print(Pw)
    H = est_homography(Pw[:,0:2], Pc)
    # print("Homography : ", H)
    z = H[2][2]
    H = H/z
    # print("Nomralised Homography : ", H)
    K_inv_H = np.matmul(np.linalg.inv(K), H)
    # print(K_inv_H)
    # Following slides: Pose from Projective Transformation
    H_dash = np.zeros((3,3))
    H_dash[:, 0] = K_inv_H[:, 0]
    H_dash[:, 1] = K_inv_H[:, 1]
    # print(K_inv_H[:, 0].shape)
    H_dash[:, 2] = np.cross(np.transpose(K_inv_H[:, 0]), np.transpose(K_inv_H[:,1]))
    # print(H_dash)
    [U, S , Vt ] = np.linalg.svd(H_dash)
    mat = np.eye(3)
    mat[2][2] = np.linalg.det(np.matmul(U, Vt))
    R = np.transpose(np.matmul(np.matmul(U, mat), Vt))
    t = -np.matmul(R, (K_inv_H[:,2]/(np.linalg.norm(K_inv_H[:,0]))))

    ##### STUDENT CODE END #####
    return R, t

# R_wc =  np.array([[ 0.876204,   -0.15612653,  0.45595071],
#  [-0.47769752, -0.40661165,  0.77876315],
#  [ 0.06380928, -0.90016191, -0.43085602]])
# t_wc =  [-0.02822914, -0.04410593, -0.71409775]
# rt = np.zeros((3,3))
# rt[:,0] = R_wc[:,0]
# rt[:,1] = R_wc[:,1]
# # rt[:,2] = R_wc[:,2]
# rt[:,2] = t_wc
# print(rt)
# K = np.array([[823.8, 0.0, 304.8],
#                   [0.0, 822.8, 236.3],
#                     [0.0, 0.0, 1.0]])
# H = np.matmul(K,rt)
# print(H)
# pw = np.array([-0.07, -0.07, 0])
# pc = np.matmul(H,pw)
# print(pc)
# print(pc[0]/pc[2])
# print(pc[1]/pc[2])                