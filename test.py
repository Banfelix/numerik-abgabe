import numpy as np
'''
N = 8
x_end = 1
x_start = -1

u = np.zeros((1, N))
u[0] = [2, 5, 4, 3, 1, 4, 2, 3]

u_left = np.zeros((1, N))
u_right = np.zeros((1, N))

dx = (x_end - x_start) / N

def minmod(a, b):
    if a * b <= 0:
        return 0
    else:
        return a if abs(a) < abs(b) else b

for i in range(1, N-1):
    a = u[0, i] - u[0, i-1]
    b = u[0, i+1] - u[0, i]
    slope = minmod(a, b)
    u_left[0, i] = u[0, i] - slope * dx / 2
    u_right[0, i] = u[0, i] + slope * dx / 2

    u_left[0, 0] = u[0, 0]
    u_right[0, 0] = u[0, 0]
    u_left[0, N-1] = u[0, N-1]
    u_right[0, N-1] = u[0, N-1]

print("dx =", dx)
print("u =", u)
print("u_left =", u_left)
print("u_right =", u_right)

w = np.zeros((N, 3))
print(w)
H_loc = np.zeros(N)
H_loc[:] = 3
print(H_loc) 
w[:, 0] = H_loc
print(w)
x = np.array([1,2,3,4,5,6,7,8])
print(x)
x[:] = x[:] * H_loc[:]
print(x)

l = x[1:-1]
print(l)

'''


for i in range(6):
    print(i)