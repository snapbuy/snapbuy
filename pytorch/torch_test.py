import torch
x = torch.rand(5, 3)
print(x)
print('Rank  of x: ', x.ndim)
print('Shape of t: ', x.shape)
print('######################################1')

import numpy as np
t = np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.], [10., 11., 12.]])
print(t)
print('Rank  of t: ', t.ndim)
print('Shape of t: ', t.shape)
print('t[0] t[1] t[-1] = ', t[0], t[1], t[-1])
print('t[2:5] t[4:-1]  = ', t[2:5], t[4:-1])
print('t[:2] t[3:]     = ', t[:2], t[3:])
print('######################################2')

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
print('######################################3')

# Braodcasting
m1 = torch.FloatTensor([[3, 3]])
m2 = torch.FloatTensor([[2, 2]])
print(m1 + m2)

m1 = torch.FloatTensor([[1, 2]])
m2 = torch.FloatTensor([3]) # [3] -> [3, 3]
print(m1 + m2)
print('######################################4')

a = torch.ones(5)
print(a)
b = a.numpy()
print(b)


a = np.ones(5)
b = torch.from_numpy(a)
np.add(a, 1, out=a)
print(a)
print(b)
print('######################################5')


# let us run this cell only if CUDA is available
# We will use ``torch.device`` objects to move tensors in and out of GPU
if torch.cuda.is_available():
    device = torch.device("cuda")          # a CUDA device object
    y = torch.ones_like(x, device=device)  # directly create a tensor on GPU
    x = x.to(device)                       # or just use strings ``.to("cuda")``
    z = x + y
    print(z)
    print(z.to("cpu", torch.double))       # ``.to`` can also change dtype together!
    print('######################################6')
