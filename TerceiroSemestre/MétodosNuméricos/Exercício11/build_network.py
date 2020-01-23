import numpy as np
from PIL import Image
import cv2

# load the image
img = cv2.imread('diagrama-servicos.jpg')

lower_gray = np.array([103, 86, 65], dtype = "uint8")
upper_gray = np.array([145, 133, 128], dtype = "uint8")

mask_gray = cv2.inRange(img, lower_gray, upper_gray)
skeleton = cv2.bitwise_and(img, img, mask=mask_gray)

cv2.imshow("images", skeleton)
cv2.waitKey(0)
cv2.destroyAllWindows()

# green nodes RGB=[30,148,88]
lower_green = np.array([80,140,20], dtype = "uint8")
upper_green = np.array([95,160,35], dtype = "uint8")

mask_green = cv2.inRange(img, lower_green, upper_green)
nodes = cv2.bitwise_and(img, img, mask=mask_green)

cv2.imshow("images", nodes)
cv2.waitKey(0)
cv2.destroyAllWindows()


# vis = np.concatenate((skeleton, nodes), axis=1)
# cv2.imwrite('out.png',vis)

h_skeleton, w_skeleton = skeleton.shape[:2]
h_nodes, w_nodes = nodes.shape[:2]

vis = np.zeros((max(h_skeleton,h_nodes), w_skeleton+w_nodes, 3), np.uint8)

vis[:h_skeleton, :w_skeleton, :3] = skeleton
vis[:h_nodes, w_skeleton:w_skeleton+w_nodes, :3] = nodes

skeleton_and_nodes = Image.fromarray(vis.astype('uint8'), mode='L')
skeleton_and_nodes.show()
