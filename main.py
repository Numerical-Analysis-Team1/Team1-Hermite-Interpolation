import sympy as sp
from sympy.utilities.lambdify import lambdify
x = sp.sympify('x')
e = sp.sympify('e')


def Lagrange_i(data_points, i):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th lagrange function
    """
    li, xi = 1, data_points[i][0]

    for j in range(len(data_points)):
        xj = data_points[j][0]

        if j != i:
            li *= (x - xj) / (xi - xj)
    return li


def Lagrange_squared(data_points, i):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th squared lagrange non-function
    """
    l2 = Lagrange_i(data_points, i) ** 2
    print("Li(x)^2 =" + str(l2))
    print()
    return l2



def Lagrange_Derivative(data_points, i):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :param i: index - int
    :return: the i-th Derivative lagrange function
    """
    dl = Lagrange_i(data_points, i).diff(x)
    print("L'(x) = " + str(dl))
    print()
    return lambdify(x, dl)


def Hermite_Interpolation(data_points):
    """
    :param data_points: list of trios (x, f(x), f'(x))
    :return: the hermite function
    """
    print("\nfor the Data points:\n" + str(data_points) + "\n")
    hermite_function = 0
    for i in range(len(data_points)):
        xi, yi, mi = data_points[i]
        la_squared = Lagrange_squared(data_points, i)
        la_derivative = Lagrange_Derivative(data_points, i)
        hermite_function += (1 - 2*(x - xi) * la_derivative(xi)) * la_squared * yi + (x - xi)*la_squared * mi
        print("H" + str(i) + "(x)=" + str(hermite_function) + "\n")

    h = hermite_function
    dh = hermite_function.diff(x)
    print("\nH(x) = " + str(h) + "\n")
    h = lambdify(x, h)
    dh = lambdify(x, dh)

    for xi, yi, mi in data_points:
        print("x=" + str(xi) + " f(x)=" + str(yi) + " f'(x)=" + str(mi) + " H(x)=" + str(h(xi)) + " H'(x)=" + str(dh(xi)))

    return hermite_function


def Drive():
    """
    getting user data and calculate H(x)
    :return: None
    """
    print("please put you Data points as follows:")
    data_points = list()
    while True:
        xi = float(input("xi: "))
        yi = float(input("yi: "))
        mi = float(input("mi: "))
        data_points.append((xi, yi, mi))

        stop = input("Continue? N - no, any other key - Yes ")
        if stop == "n" or stop == "N":
            break

    h = Hermite_Interpolation(data_points)
    print("H(x) = " + str(h))


Drive()

















