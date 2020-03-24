#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Dixon functions
"""

from mathics.builtin.base import Builtin, PostfixOperator, SympyFunction
from mathics.core.expression import Expression, Integer, Number
from mathics.core.convert import (
    sympy_symbol_prefix, SympyExpression, from_sympy)
from mathics.core.rules import Pattern
from mathics.core.numbers import dps
from mathics.builtin.scoping import dynamic_scoping
from mathics.builtin.strings import String
import sympy
import itertools
import random


##########The following functions are wrapped around their respective claases. The "apply fucntion of each class
#####calls the respective function. This is done because all the functions must belong to different classes in order to
###be called from the console. Also one class cannot be created inside an other so there is no access to functions inside a
###class from another class. This leaves us with the options to write the functions we are using more than once outside the classes

########ALL THE FUNCTIONS HAVE AS INPUTS AND OUTPUTS MATHICS OBJECTS!!!!!!!!!!!!!!!!

################################################################################
# (* The next two procedures compute the Dixon polynomial                     *)
# (* (DixonPolynomial).                                                       *)
# (* Arguments for DixonSub and DixonPolynomial consist of a list of n+1      *)
# (* polynomials in n variables (polys), a list of the n original variables   *)
# (* (xlist), and a list of the n auxiliary variables (alist).                *)
################################################################################


def DixonSub_apply(polys, xlist,alist):



    xleaf = xlist.leaves ####The .leaves option creates a python-list of mathics objects out of a mathics-list of mathics obejcts
    aleaf = alist.leaves
    polysleaf = polys.leaves.copy()####Copy is needed to ensure the first line of the matrix remains the same

    xsym = [i.to_sympy() for i in xleaf]#####.to_sympy is transforming a mathics object to a sympy object
    asym = [i.to_sympy() for i in aleaf]
    res_in = polysleaf.copy()
    res_out = [res_in.copy()]
    print('INSTALLED!!!')
    if (len(xleaf) + 1) != len(polysleaf):####Check for square matrix  same as  If[Length[dix]===Length[Transpose[dix]]
        # evaluation.message('DixonSub','not_eq')
        return String("The variable list is wrong.")
    for i in range(len(xleaf)): ###for every variable in x_list
        for j in range(len(polysleaf)):###for every polynomial
            res_in[j] = from_sympy(res_in[j].to_sympy().subs(xsym[i], asym[i])) ###subs is substituting the variables
            #####from_sympy is transforming a Sympy object to a mathics object
        res_out.append(Expression('List', *res_in))####the Expression command with the option list is buildng a mathics
    ###list out of a python iterable (list, tuple etc)
    return Expression('List', *res_out)###since the output is a 2d matrix (list within a list) we need to asseble each rows
###first and the put as elements of the outer list object




def DixonPolynomial_apply(polys, xlist,alist):

    xsym = [i.to_sympy() for i in xlist.leaves]
    asym = [i.to_sympy() for i in alist.leaves]
    res = DixonSub_apply(polys, xlist, alist)

    if not isinstance(res, Expression):####checking if the DixonSub_apply returned error
        return res
    rsym = res.leaves
    dixon_list = [j.leaves for j in rsym] ###two level Mathics list unpacked in 2 level python-list
    dixon_sym = [[j.to_sympy() for j in i] for i in dixon_list]
    dix_mat = sympy.Matrix(dixon_sym) ###2-level python list to sympy.Matrix object

    for i in range(1, len(xsym) + 1):
   #####together has poor results in sympy so simplify was used
        dix_mat[i, :] = sympy.simplify((dix_mat[i, :] - dix_mat[i - 1, :]) / (xsym[i - 1] - asym[i - 1]))
    # print(rsym[0])
    det_m = sympy.simplify(sympy.det(dix_mat))
    return from_sympy(det_m)

###(* DixonMatrix computes the classical Dixon Matrix.
def DixonMatrix_apply(polys, xlist,alist):

    res = DixonPolynomial_apply(polys, xlist, alist)
    if isinstance(res, String):
        return res

    xsym = [i.to_sympy() for i in xlist.leaves]
    asym = [i.to_sympy() for i in alist.leaves]
    xasym = asym + xsym ###concatenating the a,x lists
    dix_pol = res.to_sympy()
    if dix_pol == sympy.sympify(0):###check for zero
        return from_sympy(dix_pol)
    pol = sympy.poly(dix_pol, xasym)##create sympy.poly object
    maxdegvec = sympy.degree_list(pol, xasym)###this is MaxDegVec function in mathematica. No function was written since
    ###its built in
    # print(maxdegvec)
    res_in = []
    res_out = []
    mat_dim_col = 1
    mat_dim_row = 1
    for i in range(len(maxdegvec) // 2):
        mat_dim_row = mat_dim_row * (maxdegvec[i] + 1)###Number of rows in the matrix equal to the monomials that can be
        #formed from the a list
    for i in range(len(maxdegvec) // 2,len(maxdegvec) ):
        ###Number of collums in the matrix equal to the monomials that can be
        # formed from the x list
        mat_dim_col = mat_dim_col * (maxdegvec[i] + 1)
    A = sympy.zeros(mat_dim_row,mat_dim_col)
    #####The function lesslists is done by the command range which give a list with all the numbers less than the argument
    ####The itertools.products is essentially the cartesian product for the range lists that correspongd to different variables
    ####For example if the max degree in a,b is 1,2 we ll get two lists [0,1],[0,1,2] and the cartesian produst will be
    ###the couples [0,0],[0,1],[0,2],[1,0],[1,1],[1,2] this are the half-monomials as they represent the exponets for the
    ###a-list.
    ###The same proccess is repeated with the x list to get another half set of exponents
    for n, i in enumerate(itertools.product(*[range(i + 1) for i in maxdegvec[:len(maxdegvec) // 2]])): ##a-list exp
        for m,j in enumerate(itertools.product(*[range(i + 1) for i in maxdegvec[len(maxdegvec) // 2:]])):##x-list exp
        # res_in.append()
            monomial = sympy.sympify(1)
            for e1, x1 in zip(list(i)+list(j), xasym):###this loop builds the monomials after the two half sets are joined
                monomial = monomial * x1 ** e1
            res_out.append(pol.coeff_monomial(monomial))
            A[n,m] = pol.coeff_monomial(monomial) ##No problem with the constant term as with mathematica
    #####This if-statement is used only for a case where the matrix is a row vector or a col vector or a scalar
    ###In that case the "from_sympy" command will not return a 2d matrix but MAthematica alwyas gives 2d matrix(list of lists)
    if int(mat_dim_col)==1 | int(mat_dim_row)==1:
        res_in=from_sympy(A[0,0])
        res_out=[Expression('List', *[res_in])]
        return Expression('List', *res_out)
    return from_sympy(A)





def GaussElimination_apply(A):
    # (*Finds current pivot.First pivot in each sub-block.*)
    def find_first_pivot(A):



        ind = 0
        if A == sympy.zeros(*A.shape):
            return -1, -1 ###zero matrixs returns -1

        while A.T[ind] == sympy.sympify(0):###if we scan the matrix with one index it scans row-wise,so we eliminate the
            ###need for ZeroVectorQ function
            ind += 1
            ###the result of the integer division ind // (A.shape[0]) gives us the row of the index
            #####The remainder of the above ind % (A.shape[0]) gives us the collumn
            ###The row and collum are interchanged in the return statement because the index corresponds to the TRANSPOSE
        return [ind % (A.shape[0]), ind // (A.shape[0])]
    #(* Swaps rows Rowi and Rowj .
    def SwapRow(A, i, j):

        if i == j:
            return A
        row_i = A[i, :]
        row_j = A[j, :]
        A[j, :] = -row_i
        A[i, :] = row_j
        return A

    ##(*Replace Rowj with: Rowj + scal * Rowi.*)
    def Addrow(A, j, i, scal):

        A[j, :] = sympy.expand(A[j, :] + scal * A[i, :])
        return A
    ##(* The forward pivot.
    def pivot_for(A, i, j):

        if A[i, j] == sympy.sympify(0):
            return A
        for k in range(i + 1, A.shape[0]):
            A = Addrow(A, k, i, -A[k, j] / A[i, j])
        return A

    A_l = A.leaves
    A_ll = [j.leaves for j in A_l]
    A_s = [[j.to_sympy() for j in i] for i in A_ll]

    if len(A_s) == 1:
        return A
    Am = sympy.Matrix(A_s)

    Amt = Am ####AMt is the matrix which will be reduced in size as we procced
    top = -1 ###The variables are initialized at -1 beacuse python is ZERO-BASED and MATHEMATICA is ONE-BASED!!!!!!!
    t1 = 0
    t2 = -1
    fp = find_first_pivot(Amt)###Python cant assign in a condition as mathematica so the first pivot is calculated
    ###rigth before the while condition check. I
    # Am=SwapRow(Am, 0, 1)
    while Amt.shape[0] != 1 and fp[0] != -1:
        fp[0] = fp[0] + top + 1###+1 used because python is ZERO BASED
        fp[1] = fp[1] + t2 + 1
        top += 1
        Am1 = SwapRow(Am, top, fp[0])
        Am = Am1
        Am1 = pivot_for(Am, top, fp[1])
        Am = Am1
        Amt = Am[top + 1:, fp[1] + 1:]###no need for MinorColBlock function
        t1 = fp[0]
        t2 = fp[1]
        fp = find_first_pivot(Amt)###Python cant assign in a condition as mathematica so the first pivot is calculated
    ###rigth before the while condition check.

    return from_sympy(sympy.simplify(Am))

class DixonSub(Builtin):
    ###every class comunicates with the console with the function apply which must stated as below with the required
    ##self and evaluation objects and the input arguments

    ret_top='True'
    def apply(self, polys, xlist,alist, evaluation):
        'DixonSub[polys_, xlist_,alist_]'
        ###the line above is also needed for correct calling of the class from the console
        return DixonSub_apply(polys, xlist,alist)

class DixonPolynomial(Builtin):
    def apply(self, polys, xlist,alist, evaluation):
        'DixonPolynomial[polys_, xlist_,alist_]'
        return DixonPolynomial_apply(polys, xlist, alist)


class DixonMatrix(Builtin):
    def apply(self, polys, xlist,alist, evaluation):
        'DixonMatrix[polys_, xlist_,alist_]'
        return DixonMatrix_apply(polys, xlist,alist)


# (* ClassicalDixonResultant computes the Dixon Resultant, i.e.,             *)
# (* the determinant of the Dixon Matrix.                                    *)
#
# (* Note:  The number of polynomials, polys, MUST be one more than the      *)
# (* number of variables, vars, otherwise the determinant cannot be computed.*)

class ClassicalDixonResultant(Builtin):
    def apply(self, polys, xlist,alist, evaluation):
        'ClassicalDixonResultant[polys_, xlist_,alist_]'
        ####this is not written in a wrapper function as it is not used anywhere else

        res=DixonMatrix_apply(polys, xlist, alist)
        if isinstance(res, String):
            return res
        if res.to_sympy() == sympy.sympify(0):
            return res
        rsym = res.leaves
        dixon_list = [j.leaves for j in rsym]

        dixon_sym = [[j.to_sympy() for j in i] for i in dixon_list]
        dix_mat = sympy.Matrix(dixon_sym)
        return from_sympy(sympy.det(dix_mat))



class GaussElimination(Builtin):

    def apply(self, A, evaluation):
        'GaussElimination[A_]'
        return GaussElimination_apply(A)




class DixonResultant(Builtin):

    def apply(self, polys, xlist,alist, evaluation):
        'DixonResultant[polys_, xlist_,alist_]'
        res=DixonMatrix_apply(polys, xlist, alist)
        if isinstance(res, String):
            return res
        if res.to_sympy() == sympy.sympify(0):
            return res
        rsym = res.leaves
        dixon_list = [j.leaves for j in rsym]

        dixon_sym = [[j.to_sympy() for j in i] for i in dixon_list]
        dix_mat = sympy.Matrix(dixon_sym)
        m, n = dix_mat.shape
        rows = [i for i in range(m) if any(dix_mat[i, j] != 0 for j in range(n))]#find non zero row indices
        cols = [j for j in range(n) if any(dix_mat[i, j] != 0 for i in range(m))]#find non zero collumn indices
        dix_mat = dix_mat[rows, cols]###keep only the non zero, no need for DeleteZeroRows,DeleteZeroColumns


        res = GaussElimination_apply(from_sympy(dix_mat))###we need to transform back to the Mathics object
        ###the function returns mathics object the we need to unwrap and transform
        rsym = res.leaves
        dixon_list = [j.leaves for j in rsym]

        dixon_sym = [[j.to_sympy() for j in i] for i in dixon_list]
        dix_mat = sympy.Matrix(dixon_sym)
        m, n = dix_mat.shape
        rows = [i for i in range(m) if any(dix_mat[i, j] != 0 for j in range(n))]
        cols = [j for j in range(n) if any(dix_mat[i, j] != 0 for i in range(m))]
        dix_mat = dix_mat[rows, cols]

        ####this is the implementaion of the ProductLeadingEntries

        dix_list = []  ###since the rows will have different amount of elements we need to use nested list not matrix
        for i in range(dix_mat.shape[0]):
            dix_line=[]##the second level list
            for j in range(dix_mat.shape[1]):
                if dix_mat[i,j]!=sympy.sympify(0):###DeleteCases[m, 0, 2] we keep only the non zero elements of every line
                    dix_line.append(dix_mat[i,j])
            if len(dix_line)>0:###DeleteCases[mm, {}] we keep only the non-empty lines and append them to the outer-list
                dix_list.append(dix_line)
        prod=1
        for i in dix_list:###for every-list-line we muliply the first element
            prod=prod* i[0]
        return from_sympy(prod)



