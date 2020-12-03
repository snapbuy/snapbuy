# from https://wikidocs.net/52846, I rearranged codes to practise psersonal purpose

import torch
import numpy as np

t = np.array([[[0, 1, 2],
               [3, 4, 5]],
              [[6, 7, 8],
               [9, 10, 11]]])
ft = torch.FloatTensor(t)

print(t)
print(t.shape)
print('######################################1')
print(ft)
print(ft.shape)
print('######################################2')

torch.Size([2, 2, 3])

print(ft.view([-1, 3])) # ft라는 텐서를 (?, 3)의 크기로 변경
print(ft.view([-1, 3]).shape)

a = ([[ 0.,  1.,  2.],
        [ 3.,  4.,  5.],
        [ 6.,  7.,  8.],
        [ 9., 10., 11.]])
torch.Size([4, 3])

