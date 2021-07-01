import sympy as sp
from sympy.utilities.lambdify import lambdify
x = sp.sympify('x')


def Lagrange(points_list):
    mulls = []
    mull = 1
    for i in range(len(points_list)):
        for j in range(len(points_list)):
            if i != j:
                mull *= (x - points_list[j][0]) / (points_list[i][0] - points_list[j][0])
        mulls.append(mull)
        mull = 1
    ans = 0
    for i in range(len(points_list)):
        ans += mulls[i] * points_list[i][1]
    return ans


def hermite_interpolation(points_list):
    l = Lagrange(points_list)
    dl = l.diff()
    l = lambdify(x, l)
    dl = lambdify(x, dl)














