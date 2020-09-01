import math

def cp_x2(a, b, c):

    """ax^2+bx+c=0"""

    if a==0:

        if b==0:

            re = [0, 0] # [实数根个数, 复数根个数]

        elif type(b)==complex or type(c)==complex:

            re = [0, 1, -1*c/b]

        else:

            re = [1, 0, -1*c/b]

    else :

        det = b*b-4*a*c

        if type(det)==complex:

            xa = (-1*b+(det)**(1/2))/(2*a)

            xb = (-1*b-(det)**(1/2))/(2*a)

            re = [0, 2, xa, xb]

        else :

            if det ==0:

                re = [2, 0, -1*b/(2*a), -1*b/(2*a)]

            else :

                xa = (-1*b+(det)**(1/2))/(2*a)

                xb = (-1*b-(det)**(1/2))/(2*a)

                if det>0:

                    re = [2, 0, xa, xb]

                else :

                    re = [0, 2, xa, xb]

    return re

def cp_x3(a, b, c, d):

    """ax^3+bx^2+cx+d=0"""

    if a==0:

        re = cp_x2(b, c, d)

    else:

        A=b*b-3*a*c

        B=b*c-9*a*d

        C=c*c-3*b*d

        det = B*B-4*A*C

        if A==B==0:

            x=-1*c/b

            re = [3, 0, x, x, x]

        elif det > 0:

            ya=A*b+3*a*((-B+(det)**(1/2))/2)

            yb=A*b+3*a*((-B-(det)**(1/2))/2)

            xa=(-b-(ya**(1/3)+yb**(1/3)))/3/a

            xb=(-b+((ya**(1/3)+yb**(1/3))/2+(ya**(1/3)-yb**(1/3))/2*3**(1/2))*(-1)**(1/2))/3/a

            xc=(-b+((ya**(1/3)+yb**(1/3))/2-(ya**(1/3)-yb**(1/3))/2*3**(1/2))*(-1)**(1/2))/3/a

            re = [1, 2, xa, xb, xc]

        elif det == 0:

            K=B/A

            xa=-b/a+K

            xb=-K/2

            re = [3, 0, xa, xb, xb]

        else :

            import math

            T= (2*A*b-3*a*B)/(2*A**(3/2))

            si=math.acos(T)

            xa=(-b-2*A**(1/2)*math.cos(si/3))/(3*a)

            xb=(-b+A**(1/2)*(math.cos(si/3)+3**(1/2)*math.sin(si/3)))/(3*a)
            
            xc=(-b+A**(1/2)*(math.cos(si/3)-3**(1/2)*math.sin(si/3)))/(3*a)

            re = [3, 0, xa, xb, xc]

    return re

def cp_x4(a, b, c, d, e):
    
  """ax^4+bx^3+cx^2+dx+e=0
  费拉里法"""

    if a == 0:

        re = cp_x3(b, c, d, e)

    else :

        P=(c*c+12*a*e-3*b*d)/9

        Q=(27*a*d*d+2*c*c*c+27*b*b*e-72*a*c*e-9*b*c*d)/54

        D=(Q*Q-P*P*P)**(1/2)

        u=(Q+D)**(1/3) if abs((Q+D)**(1/3)) > abs((Q-D)**(1/3)) else (Q-D)**(1/3)

        if u==0:

            v=0

        else :

            v=P/u

        w = -0.5+3**(1/2)/2j

        mst={}

        for k in [1, 2, 3]:

            m = (b*b-8*a*c/3+4*a*(w**(k-1)*u+w**(4-k)*v))**(1/2)

            if m==0:

                s, t = b*b-8*a*c/3, 0

            else :

                s = 2*b*b-16*a*c/3-4*a*(w**(k-1)*u+w**(4-k)*v)

                t = (8*a*b*c-16*a*a*d-2*b*b*b)/m

            mst[m]=(s, t)

        m_li = mst.keys()

        m=0

        for it in m_li:

            if abs(it)>abs(m):

                m= it

        s, t = mst[m][0], mst[m][1]

        x=[]

        for i in [1, 2, 3, 4]:

            x.append((-1*b+(-1)**(math.ceil(i/2))*m+(-1)**(i+1)*(s+(-1)**(math.ceil(i/2))*t)**(1/2))/(4*a))

        re = [0, 0]

        for j in x:

            if type(j)==complex:

                if abs(j.imag) > 1e-6:

                    re[1]+=1

                    re.append(j)

                else :

                    re[0]+=1

                    re.append(j.real)

            else :

                re[0]+=1

                re.append(j)

    return re

def cp_x5(a, b, c, d, e, f):
    
    if a == 0:
        
        re = cp_x4(b, c, d,e, f)
        
    else :
        
        pass
    

def comp(li, name="x"):
    
    """规范输出"""

    re = ""

    #dictx = {1:"x1", 2:"x2", 3:"x3", 4:"x4", 5:"x5"}

    for i in range(li[0]+li[1]):

        re = re+"{}{} = {:.9f}\n".format(name, i+1, li[i+2])

    print(re)

def main():

    print("x^4-1= 0")

    x = cp_x4(1, 0, 0, 0, -1)# x^4+1=0
    
    comp(x)

    print("(x-1)(x-2)(x-3)(x-4)=0")

    y = cp_x4(1, -10, 35, -50, 24)

    comp(y)

    print("(x-1)^4=0")

    z = cp_x4_4(1, -4, 6, -4, 1)

    comp(z)


if __name__=="__main__":
    main()

