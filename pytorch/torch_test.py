import torch
x = torch.rand(5, 3)
print(x)
print('Rank  of x: ', x.ndim)
print('Shape of t: ', x.shape)
print('######################################')

import numpy as np
t = np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.], [10., 11., 12.]])
print(t)
print('Rank  of t: ', t.ndim)
print('Shape of t: ', t.shape)
print('t[0] t[1] t[-1] = ', t[0], t[1], t[-1])
print('t[2:5] t[4:-1]  = ', t[2:5], t[4:-1])
print('t[:2] t[3:]     = ', t[:2], t[3:])
print('######################################')

y = torch.FloatTensor([0., 1., 2., 3., 4., 5., 6.])
print('FloatTensor = ', y)
print(y.dim)
print(y.shape)
print(y.size)

t = torch.FloatTensor([[1., 2., 3.],
                       [4., 5., 6.],
                       [7., 8., 9.],
                       [10., 11., 12.]
                      ])
print(t)
print(t.dim())
print(t.size())
print(t[:, 1])
print(t[:, 1].size())
print(t[:, :-1])
print('######################################')

# Braodcasting
m1 = torch.FloatTensor([[3, 3]])
m2 = torch.FloatTensor([[2, 2]])
print(m1 + m2)

m1 = torch.FloatTensor([[1, 2]])
m2 = torch.FloatTensor([3]) # [3] -> [3, 3]
print(m1 + m2)
print('######################################')


a = torch.ones(5)
print(a)
b = a.numpy()
print(b)
print('######################################')
