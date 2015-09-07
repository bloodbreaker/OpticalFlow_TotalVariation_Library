import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "horn_schunck_warp_c.h":
    void HORN_SCHUNCK_MAIN ( float **f1, float **f2, float **u, float **v,
        int nx, int ny, int bx, int by, float hx,
        float hy, float m_alpha, float epsilon_d, float epsilon_s,
        float w_bright_grad,
        int num_iter_inner, int num_iter_outer,
        float n_warp_eta, float n_omega, int n_warp_levels );

    void AllocateMem2D( float ***img, int nx, int ny);

    void DisallocateMem2D( float **img, int nx, int ny);




def tv_flow_warping(np.ndarray[np.float_t, ndim=2] image1,
                    np.ndarray[np.float_t, ndim=2] image2,
                    float alpha, float w_grad_bright,float epsilon_d,
                    float epsilon_s ):



    cdef int ny = image1.shape[0]
    cdef int nx = image1.shape[1]
    cdef int bx = 2
    cdef int by = 2
    cdef int i
    cdef int j

    cdef int num_iter_inner = 1
    cdef int num_iter_outer = 15
    cdef float n_warp_eta = 0.92
    cdef float n_omega = 1.96
    cdef n_warp_levels = 200





    cdef float **f1
    cdef float **f2
    cdef float **u
    cdef float **v

    AllocateMem2D ( &f1, nx+2*bx, ny+2*by)
    AllocateMem2D ( &f2, nx+2*bx, ny+2*by)
    AllocateMem2D ( &u, nx+2*bx, ny+2*by)
    AllocateMem2D ( &v, nx+2*bx, ny+2*by)



    cdef np.ndarray[np.float_t, ndim=2] uu = np.zeros((ny,nx),dtype=np.float)
    cdef np.ndarray[np.float_t, ndim=2] vv = np.zeros((ny,nx),dtype=np.float)



    for i in range(nx):
      for j in range(ny):
        f1[i+bx][j+by] = image1[j,i]
        f2[i+bx][j+by] = image2[j,i]


    HORN_SCHUNCK_MAIN ( f1, f2, u, v, nx, ny, bx, by, 1.0, 1.0, alpha,
                        epsilon_d, epsilon_s, w_grad_bright,num_iter_inner,
                        num_iter_outer, n_warp_eta, n_omega,
                        n_warp_levels);


    for i in range(nx):
      for j in range(ny):
        uu[j,i] = u[i+bx][j+by]
        vv[j,i] = v[i+bx][j+by]


    DisallocateMem2D ( f1, nx+2*bx, ny+2*by)
    DisallocateMem2D ( f2, nx+2*bx, ny+2*by)
    DisallocateMem2D ( u, nx+2*bx, ny+2*by)
    DisallocateMem2D ( v, nx+2*bx, ny+2*by)


    print "--------------------------------------------------------------------"
    print "------------ Optical Flow Warping Scheme ---------------------------"
    print "------------ Python Module by bloodbreaker -------------------------"
    print "                                                                    "
    print "Solver Configuration:                                               "
    print "the parameter eta of warping = " + str(n_warp_eta)
    print "the maximal number of levels = " + str(n_warp_levels)
    print "the omega in the SOR solver = " + str(n_omega)
    print "number of outer loop per level = " + str(num_iter_outer)
    print "number of inner loop per level = " + str(num_iter_inner)
    print "                                                                    "
    print "Model Configuration:                                               "
    print "alpha: weight on the smoothness term = " + str(alpha)
    print "epsilon_d: contrast on the data term = " + str(epsilon_d)
    print "epsilon_s: contrast on the smoothness term = " + str(epsilon_s)
    print "w_grad_bright: weight of gradient/brightness constancy = " + str(w_grad_bright)

    print "--------------------------------------------------------------------"

    return uu,vv




