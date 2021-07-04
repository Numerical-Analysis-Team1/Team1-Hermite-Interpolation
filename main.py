import sympy as sp
from sympy.utilities.lambdify import lambdify
x = sp.sympify('x')


def Lagrange_i(data_list, i):
    """
    :param data_list: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th lagrange function
    """
    li = 1
    for j in range(len(data_list)):
        if j != i:
            li *= (x - data_list[j][0]) / (data_list[i][0] - data_list[j][0])
    return li


def Lagrange_squared(data_list, i):
    """
    :param data_list: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th squared lagrange non-function
    """
    return Lagrange_i(data_list, i) ** 2


def Lagrange_Derivative(data_list, i):
    """
    :param data_list: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th Derivative lagrange function
    """
    return lambdify(x, Lagrange_i(data_list, i).diff(x))


def Hermite_Interpolation(data_list):
    """
    :param data_list: list of trios (x, f(x), f'(x))
    :return: the hermite function
    """
    hermite_function = 0
    for i in range(len(data_list)):
        xi, yi, mi = data_list[i]
        la_squared = Lagrange_squared(data_list, i)
        la_derivative = Lagrange_Derivative(data_list, i)
        hermite_function += (1 - 2*(x - xi) * la_derivative(xi)) * la_squared * yi + (x - xi)*la_squared * mi
    return hermite_function


data = [(1, 2, 3), (2, 3, 4), (3, 4, 5)]
f = Hermite_Interpolation(data)
f = lambdify(x, f)
print(f(3))














