import math

def list(mat):

    li=[]

    for it in mat:

        for itit in it:

            li.append(itit)

    return li

def cp_x4(a, b, c, d):

    """ax^3+bx^2+cx+d=0

    return [实根个数, 复数根个数, 各个根...]

    """

    if a==0:

        if b==0:

            if c==0:

                re = [0, 0]

            else :

                re = [1, 0, -d/c]

        else:

            det = c*c-4*b*d

            if det==0:

                re = [2, 0, -c/2/b, -c/2/b]

            elif det>0:

                xa=(-c+(det)**(1/2))/2/b

                xb=(-c-(det)**(1/2))/2/b

                re = [2, 0, xa, xb]

            else:

                xa=-c/2/b+(-det)**(1/2)/2/b*(-1)**(1/2)

                xb=-c/2/b-(-det)**(1/2)/2/b*(-1)**(1/2)

                re = [0, 2, xa, xb]
                
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

class matrix:

    def __init__(self, name, row, col, list):

        self.name=name

        self.row=row

        self.col=col

        self.list=list

        self.mat=[]

        for i in range(row):

            self.mat.append(self.list[i*col:(i+1)*col])

    def eigenvalue(self):

        if self.row==1:

            re = [1, self.list[0]-1]

        elif self.row==2:

            mat=self.mat

            re = cp_x4(0, 1, -mat[0][0]-mat[1][1], mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])

        elif self.row==3:

            mat=self.mat

            a, b, c= mat[0][0], mat[0][1], mat[0][2]

            d, e, f= mat[1][0], mat[1][1], mat[1][2]

            g, h, i= mat[2][0], mat[2][1], mat[2][2]

            det = self.det()

            re = cp_x4(-1, a+e+i, c*g+h*f+b*d-a*e-a*i-e*i, det)

        elif self.row==4:

            mat=self.mat

            a, b, c, d=mat[0][0], mat[0][1], mat[0][2], mat[0][3]

            e, f, g, h=mat[1][0], mat[1][1], mat[1][2], mat[1][3]

            i, j, k, l=mat[2][0], mat[2][1], mat[2][2], mat[2][3]

            m, n, o, p=mat[3][0], mat[3][1], mat[3][2], mat[3][3]

            A=1

            B=-p-k-a-f

            C=-n*h-c*i

            D=n*h*(k+a)+c*i*(f+p)-a*f*(p+k)-p*k*(a+f)

            E=self.det()

            re = cp_x5(A, B, C, D, E)

        else:

            re = [0, 0]

        return re
    def eigenvector(self):

        eigen=self.eigenvalue()

        ei_dict=dict()

        for item in eigen[2:]:

            if type(item)!=complex:

                if item in ei_dict:

                    ei_dict[item]+=1

                else :

                    ei_dict[item]=1

        re = []

        for key in ei_dict.keys():

            temp_matrix=diag([key]*self.row)

            cop_matrix=self-temp_matrix

            coef_matrix=cop_matrix.add_col_right([0]*self.col, name="coefficient")

            so_matrix=coef_matrix.solve()

            mat = so_matrix.mat.copy()
            
            print(so_matrix)

            if self.row == 2:

                if ei_dict[key]==1:

                    re.append([-mat[0][1], mat[0][0]])

                else :

                    re.append("can't solve")

            elif self.row==3:

                if ei_dict[key]==1:

                    if mat[1][1]!=0:

                        gc=gcd(mat[0][0], mat[1][1])

                        re.append([(-1)*mat[0][2]*mat[1][1]/gc, (-1)*mat[1][2]*mat[0][0]/gc, mat[0][0]*mat[1][1]/gc])

                        print(gc)

                    elif mat[2][2]==0:

                        re.append([(-1)*mat[0][1], mat[0][0], 0])

                    else :

                        re.append("can't solve")

                elif ei_dict[key]==2:

                    re.append([(-1)*mat[0][1], mat[0][0], 0])

                    re.append([(-1)*mat[0][2], 0, mat[0][0]])

                else :

                    re.append("can't solve")
                    
            elif self.row==4:

                if ei_dict[key]==3:

                    re.append([-mat[0][1],mat[0][0], 0, 0])

                    re.append([-mat[0][2], 0,mat[0][0], 0])

                    re.append([-mat[0][3], 0, 0,mat[0][0]])

                elif ei_dict[key]==2:

                    if mat[1][1]!=0:

                        gc=gcd(mat[0][0], mat[1][1])

                        re.append([-mat[0][2]*mat[1][1]/gc, -mat[1][2]*mat[0][0]/gc,gc,0])

                        re.append([-mat[0][3]*mat[1][1]/gc, -mat[1][3]*mat[0][0]/gc,0,gc])

                    elif mat[1][2]!=0:

                        re.append([-mat[0][1], mat[0][0], 0, 0])

                        gc = gcd(mat[0][0], mat[1][2])

                        re.append([-mat[0][3]*mat[1][2]/gc, 0, -mat[1][3]*mat[0][0]/gc, 0])

                    elif mat[1][3]!=0:

                        re.append([-mat[0][1], mat[0][0], 0, 0])

                        re.append([-mat[0][2], 0, mat[0][0], 0])

                    else :
                        re.append("can't solve")

                elif ei_dict[key]==1:

                    if mat[1][1]==0:

                        re.append([-mat[0][1], mat[0][0]], 0, 0)

                    elif mat[2][2]==0:

                        gc=gcd(mat[0][0], mat[1][1])

                        re.append([-mat[0][2]*mat[1][1]/gc, -mat[1][2]*mat[0][0]/gc,gc, 0])

                    elif mat[3][3]==0:

                        gc = gcd(gcd(mat[0][0], mat[1][1]), mat[2][2])

                        re.append(-mat[0][3]*mat[1][1]*mat[2][2]/gc, -mat[1][3]*mat[0][0]*mat[2][2]/gc, -mat[2][3]*mat[0][0]*mat[1][1]/gc, gc)

                else :

                    re.append("can't solve")

            else :#self.row>4

                re.append("can't solve")

        return re
    def adjugate(self):

        pass

    def cofactor(self):

        pass

    def diagon(self):

        pass

    def rank(self):

        pass

    def add_mat(self, other):

        pass

    def add_row_below(self, in_list, ro=-2, name=""):

        mat=self.mat.copy()

        if ro =="first":

            in_ro=0

        elif ro=="end" or ro==-2:

            in_ro=self.row

        else:

            in_ro=ro

        mat.insert(in_ro, in_list)

        return matrix(name, self.row+1, self.col, list(mat))
        
    def add_col_right(self, in_list, co=-2, name=""):

        mat=self.mat.copy()

        if co =="first":

            in_co=0

        elif co=="end" or co==-2:

            in_co=self.row

        else:

            in_co=co

        for i in range(self.row):

            mat[i].insert(in_co, in_list[i])

        return matrix(name, self.row, self.col+1, list(mat))

    def __add__(self, other):

        """self+other"""

        add_list=[]

        if type(other)==matrix:

            for i in range(self.row*self.col):

                add_list.append(self.list[i]+other.list[i])

            add_mat=matrix("({}+{})".format(self.name, other.name), self.row, self.col, add_list)

        else:

            for item in self.list:

                add_list.append(item+other)

            add_mat=matrix("({}+{})".format(self.name, other), self.row, self.col, add_list)

        return add_mat
        
    def __radd__(self, other):

        """other+self"""

        add_list=[]

        if type(other)==matrix:

            for i in range(self.row*self.col):

                add_list.append(self.list[i]+other.list[i])

            add_mat=matrix("({}+{})".format(other.name, self.name), self.row, self.col, add_list)

        else:

            for item in self.list:

                add_list.append(item+other)

            add_mat=matrix("({}+{})".format(other, self.name), self.row, self.col, add_list)

        return add_mat

    def __iadd__(self, other):

        """self+=other"""

        iadd_list=[]

        for item in self.list:

            iadd_list.append(item+other)

        return matrix("({}+{})".format(self.name, other), self.row, self.col, iadd_list)
        
    def __sub__(self, other):

        """self-other"""

        sub_list=[]

        if type(other)==matrix:

            for i in range(self.row*self.col):

                sub_list.append(self.list[i]-other.list[i])

            sub_mat=matrix("({}-{})".format(self.name, other.name), self.row, self.col, sub_list)

        else:

            for item in self.list:

                sub_list.append(item-other)

            sub_mat=matrix("({}-{})".format(self.name, other), self.row, self.col, sub_list)

        return sub_mat

    def __rsub__(self, other):

        """other-self"""

        sub_list=[]

        if type(other)==matrix:

            for i in range(self.row*self.col):

                sub_list.append(self.list[i]-other.list[i])

            sub_mat=matrix("({}-{})".format(other.name, self.name), self.row, self.col, sub_list)

        else:

            for item in self.list:

                sub_list.append(other-item)

            sub_mat=matrix("({}-{})".format(other, self.name), self.row, self.col, sub_list)
            
        return sub_mat
        
    def __isub__(self, other):

        """self-=other"""

        isub_list=[]

        for item in self.list:

            isub_list.append(item-other)

        return matrix("({}-{})".format(self.name, other), self.row, self.col, isub_list)
    
    def __mul__(self, other):

        """

        func:self*other

        """

        mul_list=[]

        if type(other)==matrix:

            for i in range(self.row):

                for j in range(other.col):

                    temp=0

                    for k in range(self.col):

                        temp += self.mat[i][k]*other.mat[k][j]

                    mul_list.append(temp)

            mul_mat=matrix("({}*{})".format(self.name, other.name), self.row, other.col, mul_list)

        else :

            for item in self.list:

                mul_list.append(item*other)

            mul_mat=matrix("({}*{})".format(self.name, other), self.row, self.col, mul_list)

        return mul_mat
        
    def __rmul__(self, other):

        """

        func: other*self

        """

        mul_list=[]

        if type(other)==matrix:

            for i in range(self.row):

                for j in range(other.col):

                    temp=0

                    for k in range(self.col):

                        temp += self.mat[i][k]*other.mat[k][j]

                    mul_list.append(temp)

            mul_mat=matrix("({}*{})".format(self.name, other.name), self.row, other.col, mul_list)

        else :

            for item in self.list:

                mul_list.append(item+other)

            mul_mat=matrix("({}*{})".format(other, self.name), self.row, self.col, mul_list)

        return mul_mat
        
    def __imul__(self, other):

        """self*=other"""

        imul_list=[]

        for item in self.list:

            imul_list.append(item*other)

        return matrix("({}*{})".format(self.name, other), self.row, self.col, imul_list)

    def __floordiv__(self, divisor):

        """self//divisor"""

        floordiv_list=[]

        for item in self.list:

            floordiv_list.append(item//divisor)

            floordiv_mat=matrix("({}//{})".format(self.name, divisor), self.row, self.col, floordiv_list)

        return floordiv_mat
        
    def __truediv__(self, divisor):

        """self/divisor"""

        truediv_list=[]

        for item in self.list:

            truediv_list.append(item/divisor)

            truediv_mat=matrix("({}/{})".format(self.name, divisor), self.row, self.col, truediv_list)

        return truediv_mat

    def __div__(self, divisor):

        div_list=[]

        for item in self.list:

            div_list.append(item/divisor)

            div_mat=matrix("({}/{})".format(self.name, divisor), self.row, self.col, div_list)

        return div_mat

    def __neg__(self):

        """-self"""

        neg_list= []

        for i in range(self.row*self.col):

            neg_list.append(self.list[i]*(-1))

        neg_mat_name = self.name[1:] if self.name[0]=="-" else "-({})".format(self.name)

        neg_mat = matrix(neg_mat_name, self.row, self.col, neg_list)

        return neg_mat
        
    def __pos__(self):

        """+self"""

        return self

    def __abs__(self):

        """abs(self)"""

        abs_list= []

        for item in self.list:

            abs_list.append(item*(-1))

        abs_mat=matrix("abs({})".format(self.name), self.row, self.col, abs_list)

        return abs_mat

    def __pow__(self, other):

        """self.^other"""

        pow_list=[]

        for item in self.list:

            pow_list.append(item**other)

            pow_mat=matrix("{}.^{}".format(self.name, other), self.row, self.col, pow_list)

        return pow_mat
        
    def eye(self):

        eye_list=[]

        for i in range(self.row):

            for j in range(self.col):

                if i==j:

                    eye_list.append(1)

                else:

                    eye_list.append(0)

        return matrix("(I({}))".format(self.row), self.row, self.col, eye_list)

    def __xor__(self, other):

        """self^other"""

        if other<0:

            temp= self.inv()^(-other)

        elif other == 0:

            temp= self.eye()

        else:

            temp = self

            for i in range(1, other):

                temp *= self

            temp.name="({})^{}".format(self.name, other)

        return temp
        
    def __repr__(self):

        """print(self)"""

        re = "mat {} ({}×{}):\n".format(self.name, self.row, self.col)

        for i in range(self.row):

            for j in range(self.col):

                #re = re+str(self.list[i*self.row+j])+" "

                re = re +"{:.3f}\t".format(self.mat[i][j])

            re = re + "\n"

        return re

    def copy(self):

        list = self.list.copy()

        return matrix(self.name+"'", self.row, self.col, list)
        
    def __getitem__(self, index):

        """self[index]"""

        return self.list[index]

    def __getattr__(self, name):

        self.name=None

    def __len__(self):

        """len(self)"""

        return (row, col, row*col)

    def tran(self):

        tran_list=[]

        for i in range(self.col):

            for j in range(self.row):

                tran_list.append(self.list[i+j*self.row])

        tran_mat=matrix("tran({})".format(self.name), self.col, self.row, tran_list)

        return tran_mat
        
    def det(self):

        if self.row==1 :

            print(self.list)

            return self.list[0]

        elif self.row==2:

            return self.list[0]*self.list[3]-self.list[1]*self.list[2]

        else :

            a = self.list

            n= self.row

            sum_det = 0

            flag = 1

            b=a.copy()

            for i in range(n):

                del b[i+n:n*n:n]

                cov = b[i]

                del b[0:n:1]

                temp = matrix("", n-1, n-1, b)

                sum_det += flag*temp.det()*cov

                flag = -1*flag

                b=a.copy()

            return sum_det
            
    def del_rc(self, dr, dc):

        temp = self.list

        del temp[dr*self.row:(dr+1)*self.row:1]

        del temp[dc::self.row]

        return matrix("", self.row-1, self.col-1, temp)

    def inv(self):

        det = self.det()

        tran_inv_li=[]

        for i in range(self.row):

            for j in range(self.col):

                li=matrix("", self.row, self.col,self.list.copy())

                temp = li.del_rc(i, j)

                tran_inv_li.append(temp.det()*(-1)**(i+j)/self.det())

        inv_mat=matrix("({})^(-1)".format(self.name), self.row, self.col, tran_inv_li).tran()

        return inv_mat
        
    def change_row(self, ra, rb):

        ra-=1

        rb-=1

        ch_mat= self.mat.copy()

        ch_mat[ra], ch_mat[rb]=ch_mat[rb], ch_mat[ra]

        ch_li= list(ch_mat)

        return matrix("", self.row, self.col, ch_li)

    def multiplication(self, ro, other):

        ro-=1

        temp=[]

        mu_mat=self.mat.copy()

        for item in mu_mat[ro]:

            temp.append(item*other)

        mu_mat[ro]=temp

        mu_li=list(mu_mat)

        return matrix("", self.row, self.col, mu_li)
        
    def doubling(self, ra, rb):

        cp=self.copy()

        ra-=1

        rb-=1

        do_mat=self.mat.copy()

        index=0

        for it in do_mat[ra]:

            if it != 0:

                na=it

                break

            index+=1

        if index==self.col:

            return self

        else:

            nb=do_mat[rb][index]

            temp=[]

            for i in range(index, self.col):

                temp.append(do_mat[rb][i]*na-do_mat[ra][i]*nb)

            do_mat[rb]=temp

            if index< self.col-2:

                if temp[index+1]<0:

                    tp=[]

                    for item in temp:

                        tp.append(item*(-1))

                    do_mat[rb]=tp

            do_li=list(do_mat)

            self=cp

            return matrix("", self.row, self.col, do_li)
            
    def column(self, co):

        return self.list[co-1::self.col]

    def row(self, ro):

        return self.mat[ro-1]

    def before(self):

        be_mat=self.mat.copy()

        for i in range(self.row):

            if be_mat[i][i]==0:

                #print(be_mat)

                for j in range(i+1, self.row):

                    if be_mat[j][i]!=0:

                        be_mat[i], be_mat[j]=be_mat[j], be_mat[i]

            if be_mat[i][i]==0:

                continue

            for k in range(i+1, self.row):

                if be_mat[k][i]!=0:

                    na=be_mat[i][i]

                    nb=be_mat[k][i]

                    temp=[]

                    for l in range(self.col):

                        temp.append(be_mat[k][l]*na-be_mat[i][l]*nb)

                        #ttemp=be_mat[k][l]*na-be_mat[i][l]*nb

                        #print("{}*{}-{}*{}={}".format(be_mat[k][l], na, be_mat[i][l], nb, ttemp))

                    be_mat[k]=temp
        for n in range(self.row):

            nz=-1

            for k in range(self.col):

                if be_mat[n][k]!=0:

                    nz=k

                    break

            #print("nz", nz)

            if nz==-1 or nz==self.col-1:

                continue

            else:

                be_gcd=gcdli(be_mat[n])

                #print("gcd", be_gcd)

                tempb=[]

                flag = 1

                if be_mat[n][nz]<0:

                    flag= -1

                for it in be_mat[n]:                        tempb.append(flag*it/be_gcd)

                be_mat[n]=tempb

        be_li=list(be_mat)

        return matrix("", self.row, self.col, be_li)
        
    def after(self):

        #print("func after() begin!")

        af_mat=self.mat.copy()

        #print("dot 1")

        for i in range(self.row-1, -1, -1):

            #print("i=", i)

            nz=-1

            for j in range(self.col):#col

                #print(af_mat, i, j)

                if af_mat[i][j]!=0:

                    nz=j

                    break

            #print("nz=", nz)

            if nz==self.col-1:

                return self

            elif nz==-1:

                continue

            else:

                for k in range(nz-1, -1, -1):

                    if af_mat[k][nz]!=0:

                        na = af_mat[k][nz]

                        nb = af_mat[i][nz]

                        #print("na, nb", na, nb)

                        temp=[]

                        for l in range(self.col):

                            temp.append(af_mat[k][l]*nb-af_mat[i][l]*na)

                        #print("temp", temp)

                        ttemp=[]
                        
                        for l in range(self.col):

                            temp.append(af_mat[k][l]*nb-af_mat[i][l]*na)

                        #print("temp", temp)

                        ttemp=[]

                        tp_gcd=gcdli(temp)

                        for item in temp:

                            ttemp.append(item/tp_gcd)

                        af_mat[k]=ttemp

                #print(af_mat)

        af_li=list(af_mat)

        return matrix("", self.row, self.col, af_li)
        
    def solve(self):

        return self.before().after()

    #def ____(self):

    #    pass

def gcd(na, nb):

    na=abs(na)

    nb=abs(nb)

    while na!=0:

        if na<nb:

            na, nb=nb, na

        na=na%nb

    return nb

def gcdli(list):

    nz_li=[]

    for it in list:

        if it!=0:

            nz_li.append(abs(it))

    if nz_li==[]:

        return 1

    else:

        temp=nz_li[0]

        for i in range(1, len(nz_li)):

            temp=gcd(temp, nz_li[i])

        return temp
def change(mat, a, b):

    mat[a], mat[b]=mat[b], mat[a]

    return mat

def make_mat():

    list= []

    name= input('mat name')

    """

    row=eval(input('row:'))

    col=eval(input('col:'))

    print("details:(one by one)")

    for i in range(row*col):

        list.append(eval(input()))

    """
    
    while True:

        try:

            rc_li=input('row col').split(" ")

            row=eval(rc_li[0])

            col=eval(rc_li[1])

            if row<1 or col<1:

                print("input error, please input again")

                continue

            mat_li=input('details:').split(" ")

            for item in mat_li:

                if item!="":

                    list.append(eval(item))

            break

        except IndexError:

            print("input error, please input again")

            continue

        except NameError:

            print("input error, please input again")

            continue

    mat=matrix(name, row, col, list[:row*col])

    return mat
    
def zeros(row, col, name=""):

    if name=="":

        name="(zeros({}×{}))".format(row, col)

    zeros_list=[]

    for i in range(row*col):

        zeros_list.append(0)

    return matrix(name, row, col, zeros_list)

def ones(row, col, name=""):

    if name=="":

        name="(ones({}×{}))".format(row, col)

    ones_list=[1]*(row*col)

    #for i in range(row*col):

    #    ones_list.append(1)

    return matrix(name, row, col, ones_list)

def twos(row, col, name=""):

    if name=="":

        name="(twos({}×{})".format(row, col)

    twos_list=[2]*row*col

    #for i in range(row*col):

    #    twos_list.append(2)

    return matrix(name, row, col, twos_list)
    
def diag(list):

    dml=[]
    
    for i in range(len(list)):
    
        for j in range(len(list)):
        
            if i==j:
            
                dml.append(list[i])
                
            else:
            
                dml.append(0)
                
    return matrix("diag", len(list), len(list), dml)

def tran(a):

    tran_list=[]
    
    for i in range(a.col):
    
        for j in range(a.row):
        
            tran_list.append(a.list[i+j*a.row])
            
    tran_mat=matrix("tran({})".format(a.name), a.col, a.row, tran_list)
    
    return tran_mat
    
    #return A.tran()
    
def mat(list, name=""):

    """wrong"""

    mat_list=[]

    for it in list:

        for itit in it:

            mat_list.append(itit)

    row=len(list)

    col=len(list[0])

    return matrix(name, row, col, mat_list)

def det(A):

    n = A.row

    a = A.list

    if n==1 :

        return a[0]

    elif n==2:

        return a[0]*a[3]-a[1]*a[2]

    else :

        sum_det = 0

        flag = 1

        b=a.copy()

        for i in range(n):

            del b[i+n:n*n:n]

            cov = b[i]

            

            del b[0:n:1]

            temp = matrix("", n-1, n-1, b)

            sum_det += flag*det(temp)*cov

            flag = -1*flag

            b=a.copy()

        return sum_det

    #return A.det()
    
def inv(A):

    #inv_list=[]

    #pass

    #inv_mat=matrix("inv({})".format(A.name), A.row, A.col, inv_list)

    return A.inv()

def test():

    A=matrix('A', 3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9])

    B=matrix("B", 3, 3, [-1, 4, -2, -3, 4, 0, -3, 1, 3])

    C=matrix("C", 3, 3, [7, 4, 16, 2, 5, 8, -2, -2, -5])

    D=matrix("D", 3, 3, [2, 2, -1, 1, 3, -1, -1, -2, 2])

    E=matrix('E', 2, 2, [1, 2, 3, 4])

    F=matrix("F", 4, 4, [5, -3, 0, 9, 0, 3, 1, -2, 0, 0, 2, 0, 0, 0, 0, 2])

    print(F.eigenvalue())

    print(F.eigenvector())

if __name__=="__main__":

    test()
