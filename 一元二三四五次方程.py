import math


def cp_x2(a, b, c):
    """ax^2+bx+c=0"""
    if a == 0:
        if b == 0:
            re = [0, 0]  # [实数根个数, 复数根个数]
        elif type(b) == complex or type(c) == complex:
            re = [0, 1, -1 * c / b]
        else:
            re = [1, 0, -1 * c / b]
    else:
        det = b * b - 4 * a * c
        if type(det) == complex:
            xa = (-1 * b + det ** (1 / 2)) / (2 * a)
            xb = (-1 * b - det ** (1 / 2)) / (2 * a)
            re = [0, 2, xa, xb]
        else:
            if det == 0:
                re = [2, 0, -1 * b / (2 * a), -1 * b / (2 * a)]
            else:
                xa = (-1 * b + det ** (1 / 2)) / (2 * a)
                xb = (-1 * b - det ** (1 / 2)) / (2 * a)
                if det > 0:
                    re = [2, 0, xa, xb]
                else:
                    re = [0, 2, xa, xb]
    return re


def cp_x3(a, b, c, d):
    """ax^3+bx^2+cx+d=0"""
    if a == 0:
        re = cp_x2(b, c, d)
    else:
        A = b * b - 3 * a * c
        B = b * c - 9 * a * d
        C = c * c - 3 * b * d
        det = B * B - 4 * A * C
        if A == B == 0:
            x = -1 * c / b
            re = [3, 0, x, x, x]
        elif det > 0:
            ya = A * b + 3 * a * ((-B + det ** (1 / 2)) / 2)
            yb = A * b + 3 * a * ((-B - det ** (1 / 2)) / 2)
            xa = (-b - (ya ** (1 / 3) + yb ** (1 / 3))) / 3 / a
            xb = (-b + ((ya ** (1 / 3) + yb ** (1 / 3)) / 2 + (ya ** (1 / 3) - yb ** (1 / 3)) / 2 * 3 ** (1 / 2)) * (
                -1) ** (1 / 2)) / 3 / a
            xc = (-b + ((ya ** (1 / 3) + yb ** (1 / 3)) / 2 - (ya ** (1 / 3) - yb ** (1 / 3)) / 2 * 3 ** (1 / 2)) * (
                -1) ** (1 / 2)) / 3 / a
            re = [1, 2, xa, xb, xc]
        elif det == 0:
            K = B / A
            xa = -b / a + K
            xb = -K / 2
            re = [3, 0, xa, xb, xb]
        else:
            T = (2 * A * b - 3 * a * B) / (2 * A ** (3 / 2))
            si = math.acos(T)
            xa = (-b - 2 * A ** (1 / 2) * math.cos(si / 3)) / (3 * a)
            xb = (-b + A ** (1 / 2) * (math.cos(si / 3) + 3 ** (1 / 2) * math.sin(si / 3))) / (3 * a)
            xc = (-b + A ** (1 / 2) * (math.cos(si / 3) - 3 ** (1 / 2) * math.sin(si / 3))) / (3 * a)
            re = [3, 0, xa, xb, xc]
    return re


def cp_x4(a, b, c, d, e):
    """ax^4+bx^3+cx^2+dx+e=0
    费拉里法"""
    if a == 0:
        re = cp_x3(b, c, d, e)
    else:
        P = (c * c + 12 * a * e - 3 * b * d) / 9
        Q = (27 * a * d * d + 2 * c * c * c + 27 * b * b * e - 72 * a * c * e - 9 * b * c * d) / 54
        D = (Q * Q - P * P * P) ** (1 / 2)
        u = (Q + D) ** (1 / 3) if abs((Q + D) ** (1 / 3)) > abs((Q - D) ** (1 / 3)) else (Q - D) ** (1 / 3)
        if u == 0:
            v = 0
        else:
            v = P / u
        w = -0.5 + 3 ** (1 / 2) / 2j
        mst = {}
        for k in [1, 2, 3]:
            m = (b * b - 8 * a * c / 3 + 4 * a * (w ** (k - 1) * u + w ** (4 - k) * v)) ** (1 / 2)
            if m == 0:
                s, t = b * b - 8 * a * c / 3, 0
            else:
                s = 2 * b * b - 16 * a * c / 3 - 4 * a * (w ** (k - 1) * u + w ** (4 - k) * v)
                t = (8 * a * b * c - 16 * a * a * d - 2 * b * b * b) / m
            mst[m] = (s, t)
        m_li = mst.keys()
        m = 0
        for it in m_li:
            if abs(it) > abs(m):
                m = it
        s, t = mst[m][0], mst[m][1]
        x = []
        for i in [1, 2, 3, 4]:
            x.append((-1 * b + (-1) ** (math.ceil(i / 2)) * m + (-1) ** (i + 1) * (
                        s + (-1) ** (math.ceil(i / 2)) * t) ** (1 / 2)) / (4 * a))
        re = [0, 0]
        for j in x:
            if type(j) == complex:
                if abs(j.imag) > 1e-6:
                    re[1] += 1
                    re.append(j)
                else:
                    re[0] += 1
                    re.append(j.real)
            else:
                re[0] += 1
                re.append(j)
    return re


def cp_x5(a, b, c, d, e, f):
    """a*x^5+b*x^4+c*x^3+d*x^2+e*x+f=0
    天衍公式"""
    if a == 0:
        re = cp_x4(b, c, d, e, f)
    else:
        L = 2 * b * b - 5 * a * c
        M = 4 * b ** 3 - 15 * a * b * c + 25 * a * a * d
        N = 7 * b ** 4 + 25 * (a * c) ** 2 - 35 * a * b * b * c + 50 * a * a * b * d - 125 * e * a ** 3
        P = 4 * b ** 5 - 25 * a * c * b ** 3 + 125 * a * a * b * b * d - 625 * a * a * a * b * e + 3125 * f * a ** 4
        G = 4 * L ** 3 - 9 * M * M + 8 * L * N
        H = 10 * L * L * M - 6 * M * N + L * P
        J = 4 * L ** 4 - 4 * L * L * N + 3 * M * P
        K = M ** 4 + N ** 3 - M * N * P
        E = 2 * (G * L) ** 2 - 2 * G * G * N + 3 * G * H * M - 4 * H * H * L - G * J * L
        F = G * G * P + 3 * G * J * M - 4 * H * J * L
        A = F * F - 12 * E * E * L
        B = 6 * F ** 3 - 64 * E * E * F * L - 72 * M * E ** 3
        C = 3 * F ** 4 - 24 * E * E * F * F * L - 48 * E * E * E * F * M - 80 * E * E * E * E * L * L
        D = F * F * G + 4 * E * F * H - 4 * E * E * J
        det1 = B * B - 4 * A * C
        det2 = P * P - 4 * L ** 5
        if L == M == N == P == 0:  # (1)
            re = [5, 0, -b / (5 * a), -c / (2 * b), -d / c, -2 * e / d, -5 * f / e]
        elif L != 0 and G == H == J == 0:  # (2)
            if 7 * L * L == 4 * N:
                x1 = (-b * L - 2 * M) / (5 * a * L)
                x2 = (-2 * b * L + M) / (10 * a * L)
                re = [5, 0, x1, x2, x2, x2, x2]
            else:
                x1 = (-2 * b * L - 9 * M) / (10 * a * L)
                x3 = (-b * L + 3 * M) / (5 * a * L)
                re = [5, 0, x1, x1, x3, x3, x3]
        elif G != 0 and E == F == 0:  # (3)
            if H * H + G * J == 0:
                x1 = (-2 * b * G - 3 * H + (20 * G * G * L - 15 * H * H) ** (1 / 2)) / (10 * a * G)
                x2 = (-2 * b * G - 3 * H - (20 * G * G * L - 15 * H * H) ** (1 / 2)) / (10 * a * G)
                x3 = (-b * G + H) / (5 * a * G)
                if 20 * G * G * L - 15 * H * H > 0:
                    re = [5, 0, x1, x2, x3, x3, x3]
                else:
                    re = [3, 2, x3, x3, x3, x1, x2]
            else:
                x1 = (-b * G - H) / (5 * a * G)
                x2 = (-b * G + H + (H * H + G * J) ** (1 / 2)) / (5 * a * G)
                x4 = (-b * G + H - (H * H + G * J) ** (1 / 2)) / (5 * a * G)
                if H * H + G * J > 0:
                    re = [5, 0, x1, x2, x2, x4, x4]
                else:
                    re = [1, 4, x1, x2, x2, x4, x4]
        elif E != 0 and D == 0:  # (4)
            x4 = (-2 * b * E - F) / (10 * a * E)
            if det1 > 0:
                y1 = 10 * A * F + 15 * (-B + (B * B - 4 * A * C) ** (1 / 2)) / 2
                y2 = 10 * A * F + 15 * (-B - (B * B - 4 * A * C) ** (1 / 2)) / 2
                x1 = (-6 * b * E + 2 * F - (y1 ** (1 / 3) + y2 ** (1 / 3))) / (30 * a * E)
                x2 = (-12 * b * E + 4 * F + (y1 ** (1 / 3) + y2 ** (1 / 3))) / (60 * a * E) + (
                            3 ** (1 / 2) * (y1 ** (1 / 2) - y2 ** (1 / 2))) / (60 * a * E) * (-1)**(1/2)
                x3 = (-12 * b * E + 4 * F + (y1 ** (1 / 3) + y2 ** (1 / 3))) / (60 * a * E) - 3 ** (1 / 2) * (
                            y1 ** (1 / 2) - y2 ** (1 / 2)) / (60 * a * E) * (-1)**(1/2)
                re = [3, 2, x1, x4, x4, x2, x3]
            elif det1 < 0:
                sita = math.acos((3 * B - 4 * A * F) / (2 * A * (-5 * A) ** (1 / 2)))
                x1 = (-3 * b * E + F - (math.cos(sita / 3)) * (-5 * A) ** (1 / 2)) / (15 * a * E)
                x2 = (-6 * b * E + 2 * F + (math.cos(sita / 3) + 3 ** (1 / 2) * math.sin(sita / 3))) / (30 * a * E)
                x3 = (-6 * b * E + 2 * F + (math.cos(sita / 3) - 3 ** (1 / 2) * math.sin(sita / 3))) / (30 * a * E)
                re = [5, 0, x1, x2, x3, x4, x4]
            else:
                re = "can't solve"
        elif D != 0 and M == N == 0:  # (5)and(6)
            if det2 > 0:
                y1 = (P + (P * P - 4 * L ** 5) ** (1 / 2)) / 2
                y2 = (P - (P * P - 4 * L ** 5) ** (1 / 2)) / 2
                x1 = (-b - (y1 ** (1 / 5) + y2 ** (1 / 5))) / (5 * a)
                x2 = (-b + (1 - 5 ** (1 / 2)) * (y1 ** (1 / 5) + y2 ** (1 / 5)) / 4) / (5 * a) + (
                            10 + 2 * 5 ** (1 / 2)) * (y1 ** (1 / 5) - y2 ** (1 / 5)) / (20 * a) * (-1)**(1/2)
                x3 = (-b + (1 - 5 ** (1 / 2)) * (y1 ** (1 / 5) + y2 ** (1 / 5)) / 4) / (5 * a) - (
                            10 + 2 * 5 ** (1 / 2)) * (y1 ** (1 / 5) - y2 ** (1 / 5)) / (20 * a) * (-1)**(1/2)
                x4 = (-b + (1 + 5 ** (1 / 2)) * (y1 ** (1 / 5) + y2 ** (1 / 5)) / 4) / (5 * a) + (
                            10 - 2 * 5 ** (1 / 2)) * (y1 ** (1 / 5) - y2 ** (1 / 5)) / (20 * a) * (-1)**(1/2)
                x5 = (-b + (1 + 5 ** (1 / 2)) * (y1 ** (1 / 5) + y2 ** (1 / 5)) / 4) / (5 * a) - (
                            10 - 2 * 5 ** (1 / 2)) * (y1 ** (1 / 5) - y2 ** (1 / 5)) / (20 * a) * (-1)**(1/2)
                re = [1, 4, x1, x2, x3, x4, x5]
            elif det2 < 0:
                sita = math.acos(P / (2 * L ** (5 / 2)))
                x1 = (-b - 2 * math.cos(sita) / 5 * L ** (1 / 2)) / (5 * a)
                x2 = (-b + ((1 - 5 ** (1 / 2)) * math.cos(sita / 5) + math.sin(sita / 5) * (10 + 2 * 5 ** (1 / 2)) ** (
                            1 / 2)) / 2 * L ** (1 / 2)) / (5 * a)
                x3 = (-b + ((1 - 5 ** (1 / 2)) * math.cos(sita / 5) - math.sin(sita / 5) * (10 + 2 * 5 ** (1 / 2)) ** (
                            1 / 2)) / 2 * L ** (1 / 2)) / (5 * a)
                x4 = (-b + ((1 + 5 ** (1 / 2)) * math.cos(sita / 5) + math.sin(sita / 5) * (10 - 2 * 5 ** (1 / 2)) ** (
                            1 / 2)) / 2 * L ** (1 / 2)) / (5 * a)
                x5 = (-b + ((1 + 5 ** (1 / 2)) * math.cos(sita / 5) - math.sin(sita / 5) * (10 - 2 * 5 ** (1 / 2)) ** (
                            1 / 2)) / 2 * L ** (1 / 2)) / (5 * a)
                re = [5, 0, x1, x2, x3, x4, x5]
            else:
                re = "can't solve"
        elif D != 0 and M * N != 0 and L == K == 0:  # (7)
            z1 = (N * N / M) ** (1 / 5)
            z2 = (M * M * M / N) ** (1 / 5)
            x1 = (-b - z1 - z2) / (5 * a)
            x2 = (-b + (1 - 5 ** (1 / 2)) * z1 / 4 + (1 + 5 ** (1 / 2)) * z2 / 4) / (5 * a) + (
                        z1 * (10 + 2 * 5 ** (1 / 2)) ** (1 / 2) / 4 + z2 * (10 - 2 * 5 ** (1 / 2)) ** (1 / 2) / 4) / (
                             5 * a) * (-1)**(1/2)
            x3 = (-b + (1 - 5 ** (1 / 2)) * z1 / 4 + (1 + 5 ** (1 / 2)) * z2 / 4) / (5 * a) - (
                        z1 * (10 + 2 * 5 ** (1 / 2)) ** (1 / 2) / 4 + z2 * (10 - 2 * 5 ** (1 / 2)) ** (1 / 2) / 4) / (
                             5 * a) * (-1)**(1/2)
            x4 = (-b + (1 + 5 ** (1 / 2)) * z1 / 4 + (1 - 5 ** (1 / 2)) * z2 / 4) / (5 * a) + (
                        z1 * (10 - 2 * 5 ** (1 / 2)) ** (1 / 2) / 4 + z2 * (10 + 2 * 5 ** (1 / 2)) ** (1 / 2) / 4) / (
                             5 * a) * (-1)**(1/2)
            x5 = (-b + (1 + 5 ** (1 / 2)) * z1 / 4 + (1 - 5 ** (1 / 2)) * z2 / 4) / (5 * a) - (
                        z1 * (10 - 2 * 5 ** (1 / 2)) ** (1 / 2) / 4 + z2 * (10 + 2 * 5 ** (1 / 2)) ** (1 / 2) / 4) / (
                             5 * a) * (-1)**(1/2)
            re = [1, 4, x1, x2, x3, x4, x5]
        else:
            re = "can't solve"
    return re


def comp(li, name="x"):
    """规范输出"""
    re = ""
    if li == "can't solve":
        re = li
    else:
        for i in range(li[0] + li[1]):
            re = re + "{}{} = {:.9f}\n".format(name, i + 1, li[i + 2])
    print(re)


def main():
    print("x^4-1= 0")
    x = cp_x4(1, 0, 0, 0, -1)  # x^4+1=0
    comp(x)
    print("(x-1)(x-2)(x-3)(x-4)=0")
    y = cp_x4(1, -10, 35, -50, 24)
    comp(y)
    print("(x-1)^4=0")
    z = cp_x4(1, -4, 6, -4, 1)
    comp(z)


if __name__ == "__main__":
    main()
