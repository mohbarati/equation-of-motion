from sympy import Function, simplify, symbols

x, i, t, v, a = symbols("x i t v a")

x = Function("x")
v = Function("v")


def backward_difference(n):
    """executes the backward difference scheme

    Args:
        n (int): number of time intervals

    Returns:
        sympy exp: equation of motion
    """
    vi = [None] * (n + 1)
    xi = [None] * (n + 1)
    xi[0], vi[0] = x(0), v(0)
    for ind in range(1, n + 1):
        vi[ind] = vi[ind - 1] + a * t / n
        xi[ind] = xi[ind - 1] + vi[ind] * t / n
    return xi[-1]


def forward_difference(n):
    """executes the forward difference scheme

    Args:
        n (int): number of time intervals

    Returns:
        sympy exp: equation of motion
    """
    vi = [None] * (n + 1)
    xi = [None] * (n + 1)

    xi[0], vi[0] = x(0), v(0)
    for ind in range(0, n):

        vi[ind + 1] = vi[ind] + a * t / n
        xi[ind + 1] = xi[ind] + vi[ind] * t / n
    return xi[-1]


def equation_of_motion(n):
    """averaging the forward and backward finite difference schemes

    Args:
        n (int): number of time intervals

    Returns:
        sympy exp: more accurate representation of the final result
    """
    return simplify((forward_difference(n) + backward_difference(n)) / 2, order="none")


if __name__ == "__main__":
    n = 2
    print("number of intervals", n, ": ", equation_of_motion(n))
