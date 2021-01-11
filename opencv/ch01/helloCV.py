import sys
import cv2


print("hello", cv2.__version__)

img = cv2.imread('cat.bmp', cv2.IMREAD_COLOR) 

if img is None:
    print("Image load failed!")
    sys.exit()

cv2.namedWindow('image')
cv2.imshow('image', img)
cv2.waitKey()

cv2.destroyAllWindows()