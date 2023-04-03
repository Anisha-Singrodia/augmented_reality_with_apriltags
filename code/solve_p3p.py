import numpy as np

def Procrustes(X, Y):
    """
    Solve Procrustes: Y = RX + t

    Input:
        X: Nx3 numpy array of N points in camera coordinate (returned by your P3P)
        Y: Nx3 numpy array of N points in world coordinate
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3,) numpy array describing camera translation in the world (t_wc)
    """
    ##### CODE START #####
    A_bar = np.transpose(np.mean(Y, axis=0)) #3*1
    B_bar = np.transpose(np.mean(X, axis=0)) #3*1
    A_ = np.transpose(Y-A_bar) #3*n
    B_ = np.transpose(X-B_bar) #3*n
    ABt = np.matmul(A_, np.transpose(B_))
    [U, S , Vt ] = np.linalg.svd(ABt)
    mat = np.eye(U.shape[1])
    mat[-1][-1] = np.linalg.det(np.matmul(U, Vt))
    R = np.matmul(np.matmul(U, mat), Vt)
    t = np.transpose(A_bar) - np.matmul(R, np.transpose(B_bar))

    ##### CODE END #####
    return R, t

def P3P(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-3-Point problem, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3,) numpy array describing camera translation in the world (t_wc)

    """

    ##### CODE START #####

    # Invoke Procrustes function to find R, t
    # You may need to select the R and t that could transoform all 4 points correctly. 
    # R,t = Procrustes(Pc_3d, Pw[1:4])
    p1 = np.zeros((3,3))
    p1[:,0:-1] = Pc[1:,:]
    p1[:,2] = 1
    for  i in range(3):
        p1[i,:] = np.matmul(np.linalg.inv(K), np.transpose(p1[i,:]))

    for i in range(3):
        norm1 = np.linalg.norm(p1[i,:])
        p1[i,:] = p1[i,:]/norm1
    cos_alpha = np.dot(p1[1,:], p1[2,:])
    cos_beta = np.dot(p1[2,:], p1[0,:])
    cos_gamma = np.dot(p1[0,:], p1[1,:])

    #distances
    a = np.linalg.norm(Pw[2,:] - Pw[3,:])
    b = np.linalg.norm(Pw[3,:] - Pw[1,:])
    c = np.linalg.norm(Pw[1,:] - Pw[2,:])

    coeff = np.zeros(5)
    a2_b2 = a**2/b**2
    c2_b2 = c**2/b**2
    
    coeff[0] = (a2_b2 - c2_b2 - 1)**2 - (4*c2_b2*cos_alpha*cos_alpha)

    coeff[1] = 4*(((a2_b2-c2_b2) * (1-(a2_b2-c2_b2)) * cos_beta) - ((1-(a2_b2+c2_b2))*cos_alpha*cos_gamma) + (2*c2_b2*cos_alpha*cos_alpha*cos_beta))

    coeff[2] = 2 * ((a2_b2-c2_b2)**2 - 1 + (2 * (a2_b2 - c2_b2) * (a2_b2 - c2_b2) * cos_beta * cos_beta) + (2*(1-c2_b2)*cos_alpha*cos_alpha) - (4 * (a2_b2+c2_b2) * cos_alpha * cos_beta * cos_gamma) + (2*(1-a2_b2)*cos_gamma
    *cos_gamma))

    coeff[3] = 4*(-((a2_b2-c2_b2)*(1+(a2_b2-c2_b2))*cos_beta) + (2*a2_b2*cos_gamma*cos_gamma*cos_beta) - ((1-(a2_b2+c2_b2))*cos_alpha*cos_gamma))

    coeff[4] = ((1+a2_b2-c2_b2)**2) - (4*a2_b2*cos_gamma*cos_gamma)
    roots = np.roots(coeff)
    roots_real = roots[np.isreal(roots)].real
    roots_real = roots_real[roots_real>0]
    l2_norm = np.inf
    R = []
    t = []
    y_pred = []
    for i in range(roots_real.shape[0]):
        v = roots_real[i]
        u = ((-1+a2_b2-c2_b2)*v*v - (2*(a2_b2-c2_b2)*cos_beta*v) + 1 + a2_b2 - c2_b2)/(2*(cos_gamma - v*cos_alpha))
        s1 = np.sqrt((c*c)/(1 + (u*u) - (2*u*cos_gamma)))
        s2 = u*s1
        s3 = v*s1
        p_cam = np.zeros((3,3))
        p_cam[0,:] = p1[0,:]*s1
        p_cam[1,:] = p1[1,:]*s2
        p_cam[2,:] = p1[2,:]*s3
        r1,t1 = Procrustes(Pw[1:,:], p_cam)
        y = np.matmul(K, np.matmul(r1, np.transpose(Pw[0,:])) + t1)
        y = (y[:]/y[-1])[:-1]
        l2_norm_curr = np.linalg.norm(y - Pc[0,:])
        if l2_norm_curr < l2_norm:
            l2_norm = l2_norm_curr
            R = r1
            t = t1
            y_pred = y
    R = np.linalg.inv(R)
    t = -np.matmul(R,t)

    ##### CODE END #####

    return R, t


pc = np.array([[304.28405762, 346.36758423],[449.04196167, 308.92901611],[363.24179077, 240.77729797],[232.29425049, 266.60055542]])
pw = np.array([[-0.07, -0.07, 0.],[0.07, -0.07, 0.],[0.07, 0.07, 0.],[-0.07, 0.07, 0.]])
k = np.array([[823.8, 0., 304.8],[0., 822.8, 236.3],[0., 0., 1.]])
print(P3P(pc, pw, k))




# x = np.array([[0.12291957, 0.06189272, 0.7015967],[0.05758767, 0.00441186, 0.81126754],[-0.06603083, 0.02759465, 0.74977749]])
# y = np.array([[0.07, -0.07, 0.],[0.07, 0.07, 0.],[-0.07, 0.07, 0.]])
# print(Procrustes(x,y))
