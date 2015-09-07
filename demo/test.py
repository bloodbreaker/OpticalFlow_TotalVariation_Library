import numpy as np
import scipy.ndimage as img
from matplotlib import pyplot as plot
from matplotlib import cm
import tv_flow_warp as flow
from scipy.ndimage import filters

# read the image sequence
img1 = img.imread('f1.png').astype(np.float)
img2 = img.imread('f2.png').astype(np.float)

# check weather it is RGB
if np.size(img1.shape)==3:
  image1 = 0.299 * img1[:,:,0] + 0.587 * img1[:,:,1] + 0.114 * img1[:,:,2] 
  image2 = 0.299 * img2[:,:,0] + 0.587 * img2[:,:,1] + 0.114 * img2[:,:,2]
 
image11 = filters.gaussian_filter(image1,sigma=0.1)
image22 = filters.gaussian_filter(image2,sigma=0.1)

# compute the optical flow


# alpha -- the weight on the regularisation term;
# w_grad_bright -- weight between gradient constancy assumption and brightness
#                  constancy assumption
#epsilon_d -- the contrast parameter in the data term, which is for robust
#             statistics
#epsilon_s -- the contrast parameter in the regularisation term, which is for
#             the piece-wise smooth solution.


u,v = flow.tv_flow_warping(image11, image22, alpha=5,
                           w_grad_bright = 0.0,
                           epsilon_d = 0.01,
                           epsilon_s = 0.01)


# visualise the flow field
skipDisplay = 10 # display the flow field by every XX pixels

plot.figure()
plot.imshow(img1.astype(np.uint8),cmap = cm.Greys_r)
X,Y=np.meshgrid(np.arange(0,image1.shape[1],1),np.arange(0,image1.shape[0],1))

skip = (slice(None, None, skipDisplay), slice(None, None, skipDisplay))
norm = (u[skip]**2 + v[skip]**2+0.01)**0.5
plot.quiver(X[skip],Y[skip],u[skip],-v[skip],norm,
            pivot='tail', color='r', scale = 25,units='inches',
            headwidth=4)
plot.autoscale(tight=True)
plot.show()



