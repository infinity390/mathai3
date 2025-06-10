cmd_mode = True

import math
import copy
import itertools
from base import *
from fractions import Fraction
import parser
import math
def plog(string):
    print(string.replace("None", "not found"))
def compute(eq):
    if not eq.children:
        if eq.name == "s_e":
            return math.e
        elif eq.name == "s_pi":
            return math.pi
        return float(eq.name[2:])   
    values = [compute(child) for child in eq.children]
    ans = None
    if eq.name == "f_add":
        ans = sum(values)
    elif eq.name == "f_sub":
        ans =  values[0] - values[1]
    elif eq.name == "f_mul":
        result = 1.0
        for v in values:
            result *= v
        ans = result
    elif eq.name == "f_neg":
        ans =  -values[0]
    elif eq.name == "f_div":
        ans =  values[0] / values[1]
    elif eq.name == "f_pow":
        ans =  values[0] ** values[1]
    elif eq.name == "f_sin":
        ans =  math.sin(values[0])
    elif eq.name == "f_cos":
        ans = math.cos(values[0])
    elif eq.name == "f_tan":
        ans =  math.tan(values[0])
    elif eq.name == "f_arcsin":
        ans =  math.asin(values[0])
    elif eq.name == "f_arccos":
        ans =  math.acos(values[0])
    elif eq.name == "f_arctan":
        ans =  math.atan(values[0])
    elif eq.name == "f_log":
        ans =  math.log(values[0])
    
    if ans > 99999:
        
        x = 0/0
    
    return ans
def int_nth_root(a, n):
    if a < 0 and n % 2 == 0:
        return None
    x = round(math.pow(a, 1 / n))
    if x**n == a:
        return x
    return None
def perfect_square_root(n, nn=2):
    root = None
    if n < 0:
        root = int_nth_root(-n,nn)
        if root is not None:
            root = -root
    else:
        root = int_nth_root(n,nn)
    if root is None:
        return None
    return root if root ** nn == abs(n) else None
def pf(n, nn=2):
    if not isinstance(n, Fraction):
        n = Fraction(n)
    a = n.numerator
    b = n.denominator
    if a < 0:
        return None
    if perfect_square_root(a, nn) is not None and perfect_square_root(b, nn) is not None:
        return Fraction(perfect_square_root(a, nn),perfect_square_root(b, nn))
    return None
def int2(string):
    tmp = Fraction(string)
    if tmp.denominator == 1:
        return tmp.numerator
    return tmp
def flatten_tree(node):
    if not node.children:
        return node
    if node.name in ("f_add", "f_mul", "f_and", "f_or"):
        merged_children = []
        for child in node.children:
            flattened_child = flatten_tree(child)
            if flattened_child.name == node.name:
                merged_children.extend(flattened_child.children)
            else:
                merged_children.append(flattened_child)
        return TreeNode(node.name, merged_children)
    else:
        node.children = [flatten_tree(child) for child in node.children]
        return node
def convert_sub2neg(eq):
    def a1(eq):
        if eq.name == "f_sub":
            return TreeNode("f_add", [eq.children[0], TreeNode("f_neg", [eq.children[1]])])
        elif eq.name == "f_div":
            return TreeNode("f_mul", [eq.children[0], TreeNode("f_inv", [eq.children[1]])])
        term = TreeNode(eq.name, [])
        for child in eq.children:
            term.children.append(a1(child))
        return term
    def a2(eq):
        if eq.name == "f_neg":
            return TreeNode("f_mul", [tree_form("d_-1"), eq.children[0]])
        elif eq.name == "f_inv":
            return TreeNode("f_pow", [eq.children[0], tree_form("d_-1")])
        term = TreeNode(eq.name, [])
        for child in eq.children:
            term.children.append(a2(child))
        return term
    while "f_sub" in str_form(eq) or "f_div" in str_form(eq):
        eq = a1(eq)
    eq = flatten_tree(eq)
    while "f_neg" in str_form(eq) or "f_inv" in str_form(eq):
        eq = a2(eq) 
    return eq
def is_str_n(s, eq=None):
    s = tree_form(s)
    if s.name[:2] == "d_":
        if eq is None:
            return int2(s.name[2:])
        else:
            return int2(s.name[2:])==eq
    if eq is None:
        return None
    return False
def calc(eq):
    if eq.name[:2] == "d_":
        return Fraction(eq.name[2:])
    elif eq.name == "f_pow":
        a = calc(eq.children[0])
        b = calc(eq.children[1])
        if a is None or b is None:
            return None
        if b.denominator != 1:
            tmp2 = pf(a, b.denominator)
            if tmp2 is None:
                return None
            a = tmp2
        try:
            a = a ** b.numerator
        except:
            return None
        return a
    elif eq.name == "f_mul":
        p = 1
        for child in eq.children:
            tmp = calc(child)
            if tmp is not None:
                p *= tmp
            else:
                return None
        return p
    elif eq.name == "f_add":
        p = 0
        for child in eq.children:
            tmp  = calc(child)
            if tmp is not None:
                p += tmp
            else:
                return None
        return p
    return None
def simple(eq):
    if "f_list" in str_form(eq):
        return TreeNode(eq.name, [simple(copy.deepcopy(child)) for child in eq.children])
    if (eq.name[:2] == "f_" and eq.name != "f_add") or eq.name[:2] in ["d_","v_","s_"]:
        eq = copy.deepcopy(TreeNode("f_add", [eq]))
    if eq.name == "f_add":
        dic = {}
        num3 = Fraction(0)
        for i in range(len(eq.children)-1,-1,-1):
            tmp = calc(eq.children[i])
            if tmp is None:
                continue
            num3 += tmp
            eq.children.pop(i)
        for child in eq.children:
            num = 1
            if child.name == "f_mul":
                for i in range(len(child.children)-1,-1,-1):
                    tmp = calc(child.children[i])
                    if tmp is None:
                        continue
                    num *= tmp
                    child.children.pop(i)
                dic2 = {}
                for i in range(len(child.children)-1,-1,-1):
                    num2 = 1
                    child2 = str_form(child.children[i])
                    if child.children[i].name == "f_pow" and child.children[i].children[1].name[:2]== "d_":
                        num2 = Fraction(child.children[i].children[1].name[2:])
                        child2 = str_form(child.children[i].children[0])
                    if child2 not in dic2.keys():
                        dic2[child2] = num2
                    else:
                        dic2[child2] += num2
                newchild = TreeNode("f_mul", [])
                for key in sorted(dic2.keys()):
                    if dic2[key] == 0:
                        continue
                    if dic2[key] == 1:
                        newchild.children.append(tree_form(key))
                    else:
                        newchild.children.append(TreeNode("f_pow", [tree_form(key), tree_form("d_"+str(dic2[key]))]))
                if len(newchild.children)==1:
                    newchild = newchild.children[0]
                child = newchild
            if child.name in {"f_add", "f_mul"}:
                child = TreeNode(child.name, [tree_form(x) for x in sorted([str_form(x) for x in child.children])])
            if child.name == "f_mul" and child.children == []:
                num3 += num
            else:
                child = str_form(child)
                if child not in dic.keys():
                    dic[child] = num
                else:
                    dic[child] += num
        newchild = TreeNode("f_add", [])
        for key in sorted(dic.keys()):
            if dic[key] == 0:
                continue
            if dic[key] == 1:
                newchild.children.append(tree_form(key))
            else:
                newchild.children.append(TreeNode("f_mul", [tree_form("d_"+str(dic[key])), tree_form(key)]))
        for i in range(len(newchild.children)-1,-1,-1):
            if newchild.children[i].name[:2] == "d_":
                num3 += Fraction(newchild.children[i].name[2:])
                newchild.children.pop(i)
        if num3 != 0:
            x = tree_form("d_"+str(num3))
            newchild.children.append(x)
        if len(newchild.children)==1:
            newchild = newchild.children[0]
        elif len(newchild.children)==0:
            return tree_form("d_0")
        new = newchild
        if new.name in {"f_add", "f_mul"}:
            new = TreeNode(new.name, [tree_form(x) for x in sorted([str_form(x) for x in new.children])])
        new = flatten_tree(new)
        for i in range(len(new.children)):
            new.children[i] = simple(copy.deepcopy(new.children[i]))
        return new
def frac3(n):
    if n.denominator == 1:
        return tree_form("d_"+str(n.numerator))
    elif abs(n.numerator) == 1:
        b = tree_form("d_"+str(n.denominator*n.numerator))
        return TreeNode("f_pow", [b, tree_form("d_-1")])
    else:
        a = tree_form("d_"+str(n.numerator))
        b = tree_form("d_"+str(n.denominator))
        b = TreeNode("f_pow", [b, tree_form("d_-1")])
        return TreeNode("f_mul", [a,b])
def frac2(eq):
    if eq.name[:2] == "d_":
        return frac3(Fraction(eq.name[2:]))
    return TreeNode(eq.name, [frac2(child) for child in eq.children])
def expand_eq(eq):
    eq = simple(eq)
    eq = simple(eq)
    if "/" in str_form(eq):
        eq = frac2(eq)
        eq = simplify(eq)
        eq = simplify(eq)
    return eq
def common(eq):
    if eq.name == "f_add":
        con = []
        for child in eq.children:
            if child.name == "f_pow" and child.children[1].name[:2] == "d_" and int2(child.children[1].name[2:])<0:
                den = []
                n = int2(child.children[1].name[2:])
                if n == -1:
                    den.append(child.children[0])
                else:
                    den.append(TreeNode("f_pow", [child.children[0], tree_form("d_"+str(-n))]))
                con.append([[], den])
            elif child.name == "f_mul":
                num = []
                den = []
                for child2 in child.children:
                    if child2.name == "f_pow" and child2.children[1].name[:2] == "d_" and int2(child2.children[1].name[2:])<0:
                        n = int2(child2.children[1].name[2:])
                        if n == -1:
                            den.append(child2.children[0])
                        else:
                            den.append(TreeNode("f_pow", [child2.children[0], tree_form("d_"+str(-n))]))
                    else:
                        num.append(child2)
                con.append([num, den])
            else:
                con.append([[child], []])
        if len(con)>1 and any(x[1] != [] for x in con):
            a = TreeNode("f_add", [])
            for i in range(len(con)):
                b = TreeNode("f_mul", [])
                if con[i][0] != []:
                    b.children += con[i][0]
                for j in range(len(con)):
                    if i ==j:
                        continue
                    b.children +=  con[j][1]
                if len(b.children) == 1:
                    a.children.append(b.children[0])
                elif len(b.children) > 1:
                    a.children.append(b)
                else:
                    a.children.append(tree_form("d_1"))
            c = TreeNode("f_mul", [])
            for i in range(len(con)):
                c.children += con[i][1]
            if len(c.children)==1:
                c = c.children[0]
            c = TreeNode("f_pow", [c, tree_form("d_-1")])
            return TreeNode("f_mul", [expand2(a),c])
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(common(child))
    return arr
def group_const(eq, varname):
    def vlist(eq):
        out = []
        if eq.name[:2] == "v_":
            out.append(eq.name)
        for child in eq.children:
            out += vlist(child)
        return sorted(list(set(out)), key=lambda x: int(x[2:]))
    if eq.name =="f_mul":
        const = []
        for i in range(len(eq.children)-1,-1,-1):
            if eq.children[i].name[:2] == "d_" or varname not in vlist(eq.children[i]):
                const.append(eq.children.pop(i))
        eq.children.append(solve(product(const)))
        if len(eq.children) == 1:
            return eq.children[0]
        return eq
    if len(eq.children)==0:
        return eq
    return TreeNode(eq.name, [group_const(child, varname) for child in eq.children])
def structure(eq, s, varlist={}, varname=None):
    varlist = varlist.copy()
    if varname is not None:
        eq = group_const(copy.deepcopy(eq), varname)
    def commute(eq):
        lst = []
        eq = flatten_tree(eq)
        def start(eq):
            nonlocal lst
            if not eq.children:
                return
            if eq.name in ["f_add", "f_mul"]:
                lst.append(len(eq.children))
            for child in eq.children:
                start(child)
        item2 = None
        def start2(eq):
            nonlocal lst
            nonlocal item2
            if not eq.children:
                return eq
            coll = TreeNode(eq.name, [])
            if eq.name in ["f_add", "f_mul"]:
                for child in list(itertools.permutations(eq.children))[item2.pop(0)]:
                    coll.children.append(start2(child))
            else:
                for child in eq.children:
                    coll.children.append(start2(child))
            return coll
        out = []
        start(eq)
        for item in itertools.product(*[list(range(math.factorial(x))) for x in lst]):
            item2 = list(item).copy()
            out.append(str_form(start2(copy.deepcopy(eq))))
        return list(set(out))
    def structure2(eq, s):
        nonlocal varlist
        nonlocal varname
        def vlist(eq):
            out = []
            if eq.name[:2] == "v_":
                out.append(eq.name)
            for child in eq.children:
                out += vlist(child)
            return sorted(list(set(out)), key=lambda x: int(x[2:]))
        if s.name[:2] in ["u_", "p_"]:
            if varname is None:
                if "v_" in str_form(eq) and s.name[:2] == "p_":
                    return False
            else:
                if varname in vlist(eq) and s.name[:2] == "p_":
                    return False
            if s.name not in varlist.keys():
                varlist[s.name] = str_form(eq)
                return True
            elif varlist[s.name] == str_form(eq):
                return True
            return False
        if eq.name != s.name:
            return False
        if eq.children and s.children and len(eq.children) != len(s.children):
            return False
        for i in range(len(eq.children)):
            if not structure2(eq.children[i], s.children[i]):
                return False
        return True
    lst = commute(s)
    for item in lst:
        varlist = copy.deepcopy({})
        if structure2(copy.deepcopy(eq), tree_form(item)):
            return varlist
    return None
def common2(eq):
    eq = common(eq)
    return TreeNode(eq.name, [common2(child) for child in eq.children])
def powermerge(eq):
    def base(eq, index):
        return eq.children[index].children[0] if eq.children[index].name == "f_pow" else eq.children[index]
    def expo(eq, index):
        return eq.children[index].children[1] if eq.children[index].name == "f_pow" else tree_form("d_1")
    if eq.name == "f_mul":
        change = True
        while change:
            change = False
            l = len(eq.children)
            for i in range(2, l + 1):
                for item in itertools.combinations(range(l), i):
                    if all(base(eq, item[0]) == base(eq, item2) for item2 in item):
                        merged_exponent = TreeNode("f_add", [expo(eq, item2) for item2 in item])
                        merged_term = base(eq, item[0]) ** merged_exponent
                        indices_to_remove = sorted(item, reverse=True)
                        for j in indices_to_remove:
                            eq.children.pop(j)
                        eq.children.append(merged_term)
                        if len(eq.children) == 1:
                            eq = copy.deepcopy(eq.children[0])
                        change = True
                        break
                if change:
                    break
    return TreeNode(eq.name, [powermerge(child) for child in eq.children])
def inversehandle(eq):
    if not hasattr(eq, "children") or not eq.children:
        return eq
    tmp = structure(copy.deepcopy(eq), tree_form('f_pow\n f_pow\n  u_0\n  u_1\n u_2'), {}.copy())
    if tmp is not None:
        for key in tmp.keys():
            tmp[key] = tree_form(tmp[key])
        exponent = TreeNode("f_mul", [tmp["u_1"], tmp["u_2"]])
        base = TreeNode("f_pow", [tmp["u_0"], exponent])
        return inversehandle(base)
    if eq.name == "f_pow" and eq.children[0].name == "f_mul":
        exponent = copy.deepcopy(eq.children[1])
        ans = TreeNode("f_mul", [])
        for child in eq.children[0].children:
            ans.children.append(copy.deepcopy(TreeNode("f_pow", [child, exponent])))
        return inversehandle(ans)
    if eq.name == "f_log10":
        return inversehandle(eq.children[0].fx("log")/tree_form("d_10").fx("log"))
    if eq.name == "f_log" and eq.children[0].name == "s_e":
        return tree_form("d_1")
    tmp = structure(eq, tree_form("s_e")**(tree_form("p_0")*tree_form("u_0").fx("log")))
    if tmp is not None:
        return inversehandle(tree_form(tmp["u_0"])**tree_form(tmp["p_0"]))
    tmp = structure(eq, tree_form("s_e")**(tree_form("p_0")*tree_form("p_1")*tree_form("u_0").fx("log")))
    if tmp is not None:
        return inversehandle(tree_form(tmp["u_0"])**(tree_form(tmp["p_1"])*tree_form(tmp["p_0"])))
    if eq.name == "f_log" and eq.children[0].name == "f_pow" and eq.children[0].children[0].name == "s_e":
        return inversehandle(eq.children[0].children[1])
    if eq.name == "f_pow" and eq.children[0].name == "s_e" and eq.children[1].name == "f_log":
        return inversehandle(eq.children[1].children[0])
    if (eq.name == "f_sin" and eq.children[0].name == "f_arcsin") or (eq.name == "f_cos" and eq.children[0].name == "f_arccos") or (eq.name == "f_tan" and eq.children[0].name == "f_arctan"):
        return inversehandle(eq.children[0].children[0])
    if (eq.name == "f_cos" and eq.children[0].name == "f_arcsin") or (eq.name == "f_sin" and eq.children[0].name == "f_arccos"):
        eq2 = eq.children[0].children[0]
        eq2 = (tree_form("d_1") - eq2*eq2)**(tree_form("d_1")/tree_form("d_2"))
        return inversehandle(eq2)
    if (eq.name == "f_arcsin" and eq.children[0].name == "f_sin") or (eq.name == "f_arccos" and eq.children[0].name == "f_cos") or (eq.name == "f_arctan" and eq.children[0].name == "f_tan"):
        return inversehandle(eq.children[0].children[0])
    return TreeNode(eq.name, [flatten_tree(inversehandle(child)) for child in eq.children])
def simplify(eq, lite=False):
    def iscont(a, b, c=None):
        if c is None:
            a, b = [tree_form(x).name for x in [a,b]]
        else:
            a, b, c = [tree_form(x).name for x in [a,b,c]]
        if (c is None and all(x[:2] == "d_" for x in [a,b])) or all(x[:2] == "d_" for x in [a,b,c]):
            if c is None:
                a, b = [int(x[2:]) for x in [a,b]]
                if a == 0:
                    return None
                a = Fraction(a,1)**Fraction(b,1)
                return frac2(tree_form("d_" + str(a)))
            else:
                a, b, c = [int(x[2:]) for x in [a,b,c]]
            s = "d_"
            if a < 0 and c == -1:
                return None
            if c == -1:
                a = pf(a, abs(b))
                if a is None:
                    return None
                if b < 0:
                    a = Fraction(1,1)/a
                return frac2(tree_form("d_" + str(a)))
        return None
    if eq is None:
        return None
    if not lite:
        tmp = structure(eq, tree_form('f_pow\n u_0\n f_pow\n  u_1\n  u_2'))
        if tmp is not None:
            tmp2 = iscont(*[tmp[x] for x in sorted(tmp.keys())])
            if tmp2 is not None:
                return tmp2
        tmp = structure(eq, tree_form('f_pow\n u_0\n u_1'))
        if tmp is not None:
            tmp2 = iscont(*[tmp[x] for x in sorted(tmp.keys())])
            if tmp2 is not None:
                return tmp2
    tmp = structure(eq, tree_form('f_mul\n d_1\n u_0'))
    if tmp is not None:
        return simplify(tree_form(tmp["u_0"]), lite)
    tmp = structure(eq, tree_form('f_add\n d_0\n u_0'))
    if tmp is not None:
        return simplify(tree_form(tmp["u_0"]), lite)
    tmp = structure(eq, tree_form('f_mul\n d_0\n u_0'))
    if tmp is not None:
        return tree_form("d_0")
    tmp = structure(eq, tree_form('f_pow\n u_0\n d_0'))
    if tmp is not None:
        return tree_form("d_1")
    tmp = structure(eq, tree_form('f_pow\n u_0\n d_1'))
    if tmp is not None:
        return tree_form(tmp["u_0"])
    if all(x.name[:2] == "d_" for x in eq.children):
        if eq.name == "f_sin":
            if eq.children[0].name == "d_0":
                return tree_form("d_0")
        if eq.name == "f_cos":
            if eq.children[0].name == "d_0":
                return tree_form("d_1")
        if eq.name == "f_pow":
            if int2(eq.children[0].name[2:]) == 0 and int2(eq.children[1].name[2:]) <= 0:
                return None
            if int2(eq.children[1].name[2:]) > 0:
                n = int2(eq.children[0].name[2:]) ** int2(eq.children[1].name[2:])
                return tree_form("d_"+str(n))
            if int2(eq.children[1].name[2:]) == 0:
                return None
            if abs(int2(eq.children[0].name[2:])) == 1 and int2(eq.children[1].name[2:]) == -1:
                return eq.children[0]
        elif eq.name == "f_add":
            arr = [int2(x.name[2:]) for x in eq.children]
            return tree_form("d_"+str(sum(arr)))
        elif eq.name == "f_mul":
            arr = [int2(x.name[2:]) for x in eq.children]
            p = 1
            for item in arr:
                p = p * item
            return tree_form("d_"+str(p))
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        child = simplify(child, lite)
        if child is None:
            return None
        arr.children.append(child)
    return arr
def dowhile(eq, fx, p=False):
    while True:
        orig = copy.deepcopy(eq)
        eq = copy.deepcopy(fx(eq))
        if eq is None:
            return None
        if eq == orig:
            return orig
        elif p:
            print(eq)
def solve2(eq):
    eq = copy.deepcopy(eq)
    for item in [convert_sub2neg, expand_eq, flatten_tree, inversehandle, powermerge, simplify]:
        eq = dowhile(eq, item)
    return eq
def solve3(eq):
    eq = copy.deepcopy(eq)
    for item in [convert_sub2neg, expand_eq, flatten_tree, inversehandle, epowersplit, simplify]:
        eq = dowhile(eq, item)
    return eq
def solve(eq):
    return dowhile(eq, solve2)
def solve4(eq):
    return dowhile(eq, solve3)

def factorgen(eq):
    output = []
    if eq.name != "f_mul":
        eq = TreeNode("f_mul", [eq])
    if eq.name == "f_mul":
        for child in eq.children:
            if child.name == "f_pow":
                if child.children[0].name == "s_e":
                    output.append(child)
                    continue
                if child.children[1].name[:2] != "d_":
                    output.append(child)
                    continue
                try:
                    n = int2(child.children[1].name[2:])
                    if n < 0:
                        for i in range(-n):
                            output.append(tree_form("d_1")/child.children[0])
                    else:
                        for i in range(n):
                            output.append(child.children[0])
                except:
                    output.append(child)
            else:
                output.append(child)
    return [copy.deepcopy(simplifylite(x)) for x in output]
def substitute_val(eq, val, var="v_0"):
    eq = replace(eq, tree_form(var), tree_form("d_"+str(val)))
    return eq
def product_to_sum(eq):
    lst = factorgen(eq)
    if len(lst) == 1:
        return lst[0]
    if len(lst) == 2:
        a, b = lst
        if a.name == "f_sin" and b.name == "f_sin":
            return ((a.children[0] - b.children[0]).fx("cos") - (a.children[0] + b.children[0]).fx("cos")) / tree_form("d_2")
        elif a.name == "f_cos" and b.name == "f_cos":
            return ((a.children[0] - b.children[0]).fx("cos") + (a.children[0] + b.children[0]).fx("cos")) / tree_form("d_2")
        elif a.name == "f_sin" and b.name == "f_cos":
            return ((a.children[0] + b.children[0]).fx("sin") + (a.children[0] - b.children[0]).fx("sin")) / tree_form("d_2")
        elif a.name == "f_cos" and b.name == "f_sin":
            return ((a.children[0] + b.children[0]).fx("sin") - (a.children[0] - b.children[0]).fx("sin")) / tree_form("d_2")
    first, rest = lst[0], lst[1:]
    s = tree_form("d_0")
    eq = solve(expand2(first * product_to_sum(solve(TreeNode("f_mul", rest)))))
    if eq.name == "f_add":
        for child in eq.children:
            s += product_to_sum(child)
            s = solve(s)
    else:
        s = eq
    return s
def simplifylite(eq):
    return simplify(eq, True)
def decompose2(eq):
    lst = factorgen(eq)
    term = [x for x in lst if x.name == "f_pow" and x.children[1].name[:2]=="d_" and int(x.children[1].name[2:])<0]
    term2 = [x for x in lst if x.name != "f_pow" or x.children[1].name[:2]!="d_" or int(x.children[1].name[2:])>=0]
    lst = copy.deepcopy(term)
    term = copy.deepcopy(term2)
    
    if lst == []:
        return eq
    if term == []:
        term = tree_form("d_1")
    else:
        term = solve(expand2(product(term)))

    if not all(x.name == "f_pow" and x.children[1].name[:2] == "d_" and int(x.children[1].name[2:]) < 0 for x in lst):
        return eq
    
    v = varlist(eq)
    if len(v) != 1:
        return eq
    v = v[0]
    
    if degree(term, v)==4:
        return eq
    
    if not all(degree(x.children[0], v) <=3 for x in lst):
        return eq
    
    s = []
    final = []
    lst2 = Counter(lst)
    for item in lst2.keys():
        n = -lst2[item]
        item = copy.deepcopy(item.children[0])
        for i in range(n, 0, 1):
            final.append(item ** tree_form("d_"+str(i)))    
            if i == -1:
                s.append(item)
                continue
            s.append(item ** tree_form("d_"+str(-i)))

    if len(final) == 1:
        return eq
    
    lst2 = [x.children[0] ** tree_form("d_"+str(lst2[x])) for x in lst2.keys()]
    lst2 = [x.children[0] if x.children[1].name == "d_1" else x for x in lst2]
    lst4 = product(lst2)
    alloclst = []
    for i in range(0,26):
        if "v_"+str(i) != v:
            alloclst.append("v_"+str(i))
    used = []
    def allocate():
        nonlocal alloclst
        nonlocal used
        tmp = alloclst.pop(0)
        used.append(tmp)
        return tree_form(tmp)
    lst3 = []
    
    for index, x in enumerate(s):
        
        if degree(final[index].children[0], v) == 3:
            
            a = allocate()
            b = allocate()
            c = a*tree_form(v)+b
            lst3.append(solve((a*tree_form(v)+b)*lst4/x))
        else:
            lst3.append(solve(allocate()*lst4/x))
    lst3 = summation(lst3) - term
    lst3 = solve(expand2(lst3))
    
    used2 = {}
    if lst3.name == "f_add":
        for child in lst3.children:
            tmp = [str_form(x) for x in factorgen(child)]
            tmp2= [tree_form(x) for x in tmp if not contain(tree_form(x), tree_form(v))]
            n = tmp.count(v)
            if n not in used2.keys():
                used2[n] = product(tmp2)
            else:
                used2[n] = product(tmp2) + used2[n]
    used3 = []
    for key in used2.keys():
        used3.append(TreeNode("f_eq", [solve(used2[key]), tree_form("d_0")]))
    out2 = {}
    
    out = and0(TreeNode("f_and", used3))
    if out.name == "f_and":
        for child in out.children:
            child = child.children[0]
            v2 = varlist(child)[0]
            child = solve(child - tree_form(v2))
            child = child * tree_form("d_-1")
            out2[v2] = child
    out3 = []
    for key in sorted(out2.keys(), key=lambda x: used.index(x)):
        out3.append(out2[key])
    out4 = []
    out3 = copy.deepcopy(out3)
    count = 0
    
    for i in range(len(s)):
        if degree(copy.deepcopy(final[i]), v) == 3:
            a = out3.pop(0)
            b = out3.pop(0)
            out4.append(solve(a*tree_form(v)+b))
        else:
            out4.append(out3.pop(0))
            
    out4=  summation([product(list(x)) for x in zip(final, out4)])
    return solve(out4)
def decompose(eq):
    tmp = factorgen(eq)
    tmp1= solve(product([item for item in tmp if "v_" in str_form(item)]))
    tmp2= solve(product([item for item in tmp if "v_" not in str_form(item)]))
    return solve(decompose2(tmp1)*tmp2)
def nested_trig(eq):
    if eq.name in ["f_sin", "f_cos"] and eq.children[0].name in ["f_sin", "f_cos"]:
        return False
    if not eq.children:
        return True
    return all(nested_trig(child) for child in eq.children)
def replace_eq3(eq):
    if nested_trig(eq):
        eq = product_to_sum(eq)
    eq = replace_eq2(eq)
    eq = solve(eq)
    return eq
def replace_eq2(eq):
    if eq.name=="f_tan":
        return replace_eq2(TreeNode("f_div", [TreeNode("f_sin", [copy.deepcopy(eq.children[0])]), TreeNode("f_cos", [copy.deepcopy(eq.children[0])])]))
    if eq.name == "f_sec":
        return replace_eq2(tree_form("d_1")/eq.children[0].fx("cos"))
    if eq.name == "f_cosec":
        return replace_eq2(tree_form("d_1")/eq.children[0].fx("sin"))
    if eq.name == "f_cot":
        return replace_eq2(eq.children[0].fx("cos")/eq.children[0].fx("sin"))
    def isneg(eq):
        eq = tree_form(eq)
        if eq.name[:2] != "d_":
            return False
        if int(eq.name[2:]) >= 0:
            return False
        return True
    tmp = structure(eq, (tree_form("u_0")*tree_form("p_0")).fx("sin"))
    if tmp is not None and isneg(tmp["p_0"]):
        return tree_form("d_-1")*(tree_form(tmp["u_0"])*tree_form(tmp["p_0"])*tree_form("d_-1")).fx("sin")
    tmp = structure(eq, (tree_form("u_0")*tree_form("p_0")).fx("cos"))
    if tmp is not None and isneg(tmp["p_0"]):
        return (tree_form(tmp["u_0"])*tree_form(tmp["p_0"])*tree_form("d_-1")).fx("cos")
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(replace_eq2(child))
    return arr
def rmdeg(eq):
    if eq.name == "f_rad":
        return solve(eq.children[0]*tree_form("s_pi")/tree_form("d_180"))
    return TreeNode(eq.name, [rmdeg(child) for child in eq.children])
def replace_eqint(eq):
    return product([x if x.name == "f_pow" and x.children[1].name[:2] == "d_" and int(x.children[1].name[2:]) < 0 else replace_eq(x) for x in factorgen(eq)])
def addition_trig(eq, depth):
    
    if eq.name in ["f_sin", "f_cos"] and eq.children[0].name == "f_add" and len(eq.children[0].children)==2:
        if eq.name == "f_sin":
            
            a = TreeNode("f_sin", [eq.children[0].children[0]])
            b = TreeNode("f_cos", [eq.children[0].children[1]])
            c = TreeNode("f_cos", [eq.children[0].children[0]])
            d = TreeNode("f_sin", [eq.children[0].children[1]])
            return replace_eq(a*b+c*d, depth), True
        elif eq.name == "f_cos":
            a = TreeNode("f_cos", [eq.children[0].children[0]])
            b = TreeNode("f_cos", [eq.children[0].children[1]])
            c = TreeNode("f_sin", [eq.children[0].children[0]])
            d = TreeNode("f_sin", [eq.children[0].children[1]])
            return replace_eq(a*b-c*d, depth), True
    return eq, False
def replace_eq(eq, depth=1):
    if eq.name[2:] in ["tan", "sin", "cos"]:
        eq = TreeNode(eq.name, [solve(eq.children[0])])
    tmp = structure(eq, (tree_form("u_0")*tree_form("d_-1")).fx("sin"))
    if tmp is not None:
        return tree_form("d_-1")*tree_form(tmp["u_0"]).fx("sin")
    tmp = structure(eq, (tree_form("u_0")*tree_form("d_-1")).fx("cos"))
    if tmp is not None:
        return tree_form(tmp["u_0"]).fx("cos")
    tmp = structure(eq, tree_form('f_add\n f_cos\n  u_0\n f_mul\n  d_-1\n  f_sin\n   u_0'))
    if tmp is not None:
        var = tree_form(tmp["u_0"])
        return tree_form("d_2")**(tree_form("d_2")**tree_form("d_-1"))*(var+tree_form("s_pi")/tree_form("d_4")).fx("cos")


    if eq.name == "f_arctan":
        if eq.children[0].name == "d_0":
            return tree_form("d_0")

    if eq.name == "f_log":
        if eq.children[0].name == "d_1":
            return tree_form("d_0")
    if eq.name[2:] in ["cos", "sin"]:
        eq2 = eq.children[0]
        pi2 = tree_form("s_pi")/tree_form("d_2")
        eq3 = solve(expand2(eq2/pi2))
        
        if eq3.name[:2] == "d_":
            v = int(eq3.name[2:])
            if "cos" == eq.name[2:]:
                ans = 0
                if v % 4 == 0:
                    ans = 1
                elif v%4 == 2:
                    ans = -1
                return tree_form("d_"+str(ans))
            else:
                ans = 0
                if v % 4 == 1:
                    ans = 1
                elif v%4 == 3:
                    ans = -1
                return tree_form("d_"+str(ans))
        else:
            eq3 = solve(expand2(eq3 / tree_form("d_2")))
            ans = [parser.take_input("3^(1/2)/2"), parser.take_input("1/2^(1/2)"), parser.take_input("1/2")]
            
            for index, item in enumerate([parser.take_input("1/3"), parser.take_input("1/4"), parser.take_input("1/6")]):
                if solve(expand2(item - eq3)) == tree_form("d_0"):
                    if "sin" == eq.name[2:]:
                        return ans[index]
                    elif "cos" == eq.name[2:]:
                        return ans[2-index]
            
            if depth>0 and diff(copy.deepcopy(eq3)) == tree_form("d_0"):
                c = 180*compute(eq3)
                item = None
                sub = None
                sgn = parser.take_input("1")
                if c <0:
                    sgn = parser.take_input("-1")
                    c = -c
                d= str(int(c//90))
                
                if d != "0":
                    sub = parser.take_input("pi/2")*parser.take_input(d)
                    eq3 = solve(expand2(sgn * eq3))
                    item =[solve(expand2(tree_form("s_pi")*eq3 - sub)) + sub]
                else:
                    item = []
                    eq = TreeNode(eq.name, [solve(expand2(sgn * tree_form("s_pi") * eq3))])
                    if eq.name[2:] == "sin":
                        eq = solve(expand2(sgn * eq))
                    
                
                for item2 in item:
                    eq4 = TreeNode(eq.name, [item2])
                    eq4, done = addition_trig(copy.deepcopy(eq4), 0)
                    if eq.name[2:] == "sin":
                        eq4 = solve(expand2(sgn * eq4))
                    if done:
                        return eq4
                    
    eq, _ = addition_trig(eq, depth)
    
    tmp = structure(eq, tree_form('f_sin\n f_mul\n  p_0\n  u_0'))
    if tmp is not None and tree_form(tmp["p_0"]).name[:2] == "d_" and int(tree_form(tmp["p_0"]).name[2:])%2 == 0:
        var = tree_form(tmp["u_0"])*tree_form(tmp["p_0"])/tree_form("d_2")
        return solve(tree_form("d_2")*var.fx("sin")*var.fx("cos"))
    '''
    fac = factorgen(eq)
    done = None
    for item in itertools.combinations(range(len(fac)), 2):
        if (fac[item[0]].name == "f_sin" and fac[item[1]].name == "f_cos") or\
           (fac[item[0]].name == "f_cos" and fac[item[1]].name == "f_sin"):
            if solve(expand2(fac[item[0]].children[0] - fac[item[1]].children[0])) == tree_form("d_0"):
                done = item
                break
    if done is not None:
        return solve(product([item for index, item in enumerate(fac) if index not in done]) * (fac[done[0]].children[0]*tree_form("d_2")).fx("sin")/ tree_form("d_2"))
    '''
    tmp = structure(eq, tree_form('f_cos\n f_mul\n  p_0\n  u_0'))
    if tmp is not None and tree_form(tmp["p_0"]).name[:2] == "d_" and int(tree_form(tmp["p_0"]).name[2:])%2 == 0:
        var = tree_form(tmp["u_0"])*tree_form(tmp["p_0"])/tree_form("d_2")
        return solve(var.fx("cos")*var.fx("cos")-var.fx("sin")*var.fx("sin"))
    
    arr = TreeNode(eq.name, [])
    for child in eq.children:
        arr.children.append(replace_eq(child, depth))
    return arr
def expand(eq):
    if eq.name == "f_mul":
        addchild = [[child2 for child2 in child.children] for child in eq.children if child.name == "f_add"]
        otherchild = [child for child in eq.children if child.name != "f_add"]
        if len(otherchild) == 1 and otherchild[0].name[:2] == "d_" and len(addchild) == 1 and isinstance(addchild, list):
            add= tree_form("d_0")
            for item in addchild[0]:
                add += otherchild[0]*item
            return add
    return TreeNode(eq.name, [expand(child) for child in eq.children])
def addm(m1, m2):
    m1, m2 = copy.deepcopy(m1), copy.deepcopy(m2)
    if m1.name != "f_list":
        return copy.deepcopy(TreeNode("f_add", [m1, m2])), tree_form("d_0")
    coll = TreeNode("f_list", [])
    for i in range(len(m1.children)):
        coll.children.append(addm(m1.children[i], m2.children[i])[0])
    return coll, m2
def scale(m, number):
    if m.name != "f_list":
        return m*number
    return TreeNode("f_list", copy.deepcopy([scale(child, number) for child in m.children]))
def conv_mat(eq):
    if eq.name == "f_list":
        return [conv_mat(child) for child in eq.children]
    else:
        return eq
def conv_list(eq):
    if isinstance(eq, list):
        return TreeNode("f_list", [conv_list(child) for child in eq])
    else:
        return eq
def matrix_multiply(matrix1, matrix2):
    matrix1, matrix2= copy.deepcopy(matrix1), copy.deepcopy(matrix2)
    if len(matrix1[0]) != len(matrix2):
        return None
    output = []
    for i in range(len(matrix1)):
        element = []
        for j in range(len(matrix2[0])):
            coll = TreeNode("f_add", [])
            for k in range(len(matrix2)):
                coll.children.append(TreeNode("f_mul", [matrix1[i][k], matrix2[k][j]]))
            element.append(coll)
        output.append(element)
    return output
def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]
def convert_nested_list(lst):
    if isinstance(lst, int):
        return tree_form("d_"+str(lst))
    if isinstance(lst, list):
        return [convert_nested_list(item) for item in lst]
    return lst
def determinant(matrix):
    matrix = convert_nested_list(matrix)
    return determinant2(matrix)
def determinant2(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    det = tree_form("d_0")
    for col in range(n):
        sub_matrix = [row[:col] + row[col+1:] for row in matrix[1:]]
        det += tree_form("d_"+str(-1 ** col)) * matrix[0][col] * determinant2(sub_matrix)
    return det
def circumcenter(matrix):
    x1, y1 = matrix[0]
    x2, y2 = matrix[1]
    x3, y3 = matrix[2]
    a = x1**tree_form("d_2")+y1**tree_form("d_2")
    b = x2**tree_form("d_2")+y2**tree_form("d_2")
    c = x3**tree_form("d_2")+y3**tree_form("d_2")
    n1 = determinant([[a,y1,1],[b,y2,1],[c,y3,1]])
    n2 = determinant([[x1,a,1],[x2,b,1],[x3,c,1]])
    d = determinant([[x1,y1,1],[x2,y2,1],[x3,y3,1]])
    two = tree_form("d_2")
    return [n1/(two*d),n2/(two*d)]
def midpoint(matrix):
    x1, y1 = matrix[0]
    x2, y2 = matrix[1]
    return [(x1+x2)/tree_form("d_2"),(y1+y2)/tree_form("d_2")]
def matrix_simp2(eq):
    if eq.name == "f_transpose":
        eq.children[0] = conv_mat(eq.children[0])
        eq.children[0] = transpose(eq.children[0])
        return conv_list(eq.children[0])
    if eq.name == "f_circumcenter":
        eq.children[0] = conv_mat(eq.children[0])
        eq.children[0] = circumcenter(eq.children[0])
        return conv_list(eq.children[0])
    if eq.name == "f_midpoint":
        eq.children[0] = conv_mat(eq.children[0])
        eq.children[0] = midpoint(eq.children[0])
        return conv_list(eq.children[0])
    if eq.name[:-1] == "f_point":
        eq.children[0] = conv_mat(eq.children[0])
        eq.children[0] = eq.children[0][int(eq.name[-1])-1]
        return conv_list(eq.children[0])
    if eq.name == "f_mag":
        return summation([(eq.children[0].children[1].children[i]-eq.children[0].children[0].children[i])**tree_form("d_2") for i in range(len(eq.children[0].children[0].children))])
    return TreeNode(eq.name, [matrix_simp2(child) for child in eq.children])
def matrix_simp(eq):
    eq = copy.deepcopy(eq)
    if eq.name == "f_eq":
        return TreeNode(eq.name, [matrix_simp(eq.children[0]), eq.children[1]])
    if eq.name == "f_and":
        return TreeNode(eq.name, [matrix_simp(child) for child in eq.children])
    if eq.name == "f_add":
        if eq.children[0].name != "f_list":
            eq.children[0] = matrix_simp(eq.children[0])
        collect = copy.deepcopy(eq.children[0])
        for i in range(1,len(eq.children)):
            if eq.children[i].name != "f_list":
                eq.children[i] = matrix_simp(eq.children[i])
            collect = addm(collect, eq.children[i])[0]
        return collect
    if eq.name == "f_mul" and any(child.name == "f_list" for child in eq.children):
        collect = copy.deepcopy([child for child in eq.children if child.name == "f_list"][0])
        collect = conv_mat(collect)
        index = [index for index, child in enumerate(eq.children) if child.name == "f_list"][0]
        for i in range(len(eq.children)):
            if i == index:
                continue
            if "f_list" in str_form(eq.children[i]):
                collect = matrix_multiply(collect, conv_mat(eq.children[i]))
            else:
                collect = scale(conv_list(collect), eq.children[i])
                collect = conv_mat(collect)
        return conv_list(collect)
    if eq.name == "f_pow":
        out = TreeNode("f_mul", [])
        for i in range(int(eq.children[1].name[2:])):
            out.children.append(copy.deepcopy(eq.children[0]))
        return matrix_simp(out)
    return eq

def expand2(eq):
    eq = copy.deepcopy(eq)
    if eq.name == "f_mul" or eq.name == "f_pow":
        if eq.name == "f_pow":
            eq = copy.deepcopy(TreeNode("f_pow", [eq]) )
        ac = []
        addchild = []
        for child in eq.children:
            tmp5 = factorgen(child)
            ac += tmp5
        tmp3 = []
        for child in ac:
            tmp2 = []
            if child.name == "f_add":
                if child.children != []:
                    for child2 in child.children:
                        tmp2.append(child2)
                else:
                    tmp2 = [child]
            else:
                tmp3.append(child)
            if tmp2 != []:
                addchild.append(tmp2)
        tmp4 = tree_form("d_1")
        for item in tmp3:
            tmp4 = tmp4 * item
        addchild.append([tmp4])
        def flatten(lst):
            flat_list = []
            for item in lst:
                if isinstance(item, list) and item == []:
                    continue
                if isinstance(item, list):
                    flat_list.extend(flatten(item))
                else:
                    flat_list.append(item)
            return flat_list
        if isinstance(addchild, list) and len(flatten(addchild))>0:
            add= tree_form("d_0")
            for item in itertools.product(*addchild):
                mul = tree_form("d_1")
                for item2 in item:
                    mul = mul * item2
                    mul = dowhile(mul, simplifylite)
                add = add + mul
                add = dowhile(add, simplifylite)
            eq = copy.deepcopy(add)
    return TreeNode(eq.name, [expand2(child) for child in eq.children])
def numdem(equation):
    num = tree_form("d_1")
    den = tree_form("d_1")
    for item in factorgen(equation):
        t = item
        if t.name == "f_pow" and "v_" not in str_form(t.children[1]) and "f_" not in str_form(t.children[1]) and int2(t.children[1].name[2:]) < 0:
            den = den*item
        else:
            num = num*item
    return [solve(num), solve(tree_form("d_1")/den)]
                
def diff(eq):
    eq = copy.deepcopy(eq)
    eq = dowhile(eq, simplifylite)
    eq = solve(eq)
    if "v_" not in str_form(eq):
        return tree_form("d_0")
    if eq.name == "f_add":
        add = tree_form("d_0")
        for child in eq.children:
            add += diff(child)
        return add
    elif eq.name == "f_pow" and eq.children[0].name == "s_e":
        return diff(eq.children[1])*eq
    elif eq.name == "f_tan":
        return diff(eq.children[0])/(eq.children[0].fx("cos")*eq.children[0].fx("cos"))
    elif eq.name == "f_log":
        return diff(eq.children[0])*(tree_form("d_1")/eq.children[0])
    elif eq.name == "f_arcsin":
        return diff(eq.children[0])/(tree_form("d_1")-eq.children[0]*eq.children[0])**parser.take_input("1/2")
    elif eq.name == "f_arccos":
        return tree_form("d_-1")*diff(eq.children[0])/(tree_form("d_1")-eq.children[0]*eq.children[0])**parser.take_input("1/2")
    elif eq.name == "f_arctan":
        return diff(eq.children[0])/(tree_form("d_1")+eq.children[0]*eq.children[0])
    elif eq.name == "f_pow" and "v_" in str_form(eq.children[1]):
        a, b = eq.children
        return a**b * ((b/a) * diff(a) + a.fx("log") * diff(b))
    elif eq.name == "f_mul":
        add = tree_form("d_0")
        for i in range(len(eq.children)):
            new = copy.deepcopy(eq)
            new.children.pop(i)
            if len(new.children)==1:
                new = new.children[0]
            add += diff(dowhile(eq.children[i], simplifylite))*new
        add = dowhile(add, simplifylite)
        return add
    elif eq.name == "f_sin":
        eq.name = "f_cos"
        return diff(eq.children[0])*eq
    elif eq.name == "f_cos":
        eq.name = "f_sin"
        return tree_form("d_-1")*diff(eq.children[0])*eq
    elif eq.name[:2] == "v_":
        return TreeNode("f_dif", [eq])
    elif eq.name == "f_pow" and "v_" not in str_form(eq.children[1]):
        base, power = eq.children
        dbase = diff(base)
        b1 = power - tree_form("d_1")
        bab1 = TreeNode("f_pow", [base, b1])
        return power * bab1 * dbase
def diffx2(equation, var="v_0"):
    if equation.name == "f_dif":
        if equation.children[0].name == var:
            return tree_form("d_1")
        return tree_form("d_0")
    return TreeNode(equation.name, [diffx2(child, var) for child in equation.children])
def diffany(eq):
    if eq.name == "f_dif":
        if eq.children[0].name[:2] == "v_":
            return tree_form("d_1")
    return TreeNode(eq.name, [diffany(child) for child in eq.children])
def diffx(equation, var="v_0"):
    equation = diff(equation)
    equation = diffx2(equation, var)
    return solve(equation)
def diffany2(equation):
    equation = diff(equation)
    equation = diffany(equation)
    equation = solve(equation)
    return equation

#eq = parser.take_input("log(x+1)")
#print([str_form(solve(eq))])
def approx(eq, var):
    n, d= numdem(eq)
    n, d=  replace_eq2(n), replace_eq2(d)
    n, d = solve(n), solve(d)
    n, d = expand(n), expand(d)
    out = []
    for equation in [n, d]:
        for item in factorgen(equation):
            tmp = structure(item, tree_form('f_sin\n u_0'), {},"v_0")
            if tmp is not None:
                for key in tmp.keys():
                    tmp[key] = tree_form(tmp[key])
                item2 = substitute_val(tmp["u_0"], 0, var.name)
                if tree_form("d_0") == expand3(solve(item2)):
                    equation = equation/item
                    equation = equation*tmp["u_0"]
                    break
                elif tree_form("d_0") == expand3(solve(tree_form("s_pi") - item2)):
                    equation = equation/item
                    equation = equation*(tree_form("s_pi") - tmp["u_0"])
                    break
            tmp = structure(item, tree_form('f_add\n d_-1\n f_pow\n  p_0\n  u_0'), {},"v_0")
            if tmp is not None:
                for key in tmp.keys():
                    tmp[key] = tree_form(tmp[key])
                
                item2 = substitute_val(tmp["u_0"], 0, var.name)
                item2 = expand4(solve(item2))
                if tree_form("d_0") == item2:
                    equation = equation/item
                    equation = solve(equation*tmp["u_0"]*tmp["p_0"].fx("log"))
                    break
            tmp = structure(item, tree_form('f_log\n f_add\n  d_1\n  u_0'), {},"v_0")
            if tmp is not None:
                for key in tmp.keys():
                    tmp[key] = tree_form(tmp[key])
                
                item2 = substitute_val(tmp["u_0"], 0, var.name)
                item2 = expand4(solve(item2))
                if tree_form("d_0") == item2:
                    equation = equation/item
                    equation = solve(equation*tmp["u_0"])
                    break
            tmp = structure(item, tree_form('f_cos\n u_0'), {},"v_0")
            if tmp is not None:
                for key in tmp.keys():
                    tmp[key] = tree_form(tmp[key])
                item2 = substitute_val(item, 0, var.name)
                
                if tree_form("d_0") == expand3(solve(item2)):
                    
                    equation = equation/item
                    equation = equation*(tree_form("d_1") - tmp["u_0"]**tree_form("d_2"))
                    break
                
        equation = solve(equation)
        out.append(equation)
    return solve(out[0]/out[1])
def subslimit(equation, var):
    equation2 = replace(equation, var, tree_form("d_0"))
    try:
        compute(equation2)
        return expand4(solve(equation2))
    except:
        return None
def lhospital(equation, var):
    
    e = replace(equation, var, tree_form("d_0"))
    try:
        compute(e)
    except:
        e = None
    if e is None:
        n, d = numdem(equation)
        
        ans1 = subslimit(expand4(solve(n)), var)
        ans2 = subslimit(expand4(solve(d)), var)
        
        if ans1 is not None and ans2 is not None and ans1 == tree_form("d_0") and ans2 == tree_form("d_0"):
            
            g, h = diffx(solve(n), var.name), diffx(solve(d), var.name)
            
            equation2 = copy.deepcopy(expand4(solve(g))/expand4(solve(h)))
            equation2 = solve(equation2)
            return equation2
        return None
    return None

def collec(eq, var):
    output = []
    def collect2(eq):
        if eq.name in ["f_pow", "f_sin", "f_cos", "f_arcsin"] and eq.children[0].name[:2] != "v_" and var in str_form(eq.children[0]):
            output.append(str_form(eq.children[0]))
        if eq.name == "f_pow" and eq.children[0].name == "s_e" and "v_" in str_form(eq):
            if eq.children[1].name[:2] != "v_":
                output.append(str_form(eq.children[1]))
            output.append(str_form(eq))
        for child in eq.children:
            collect2(child)
    def collect3(eq):
        if eq.name in ["f_sin", "f_cos"]:
            output.append(str_form(eq.children[0].fx("cos")))
        for child in eq.children:
            collect3(child)  
    collect2(eq)
    if output == []:
        collect3(eq)
    tmp = sorted(output, key=lambda x: len(x))
    tmp = [solve(tree_form(x)) for x in tmp]
    
    return tmp

constant_term= 0
def ccc():
    global constant_term
    constant_term += 1
    return tree_form("v_c"+str(constant_term-1))
def byparts(f, g, depth, var, data, bypart, sp, su):
    eq = integratex(g, depth-1, var, data+[TreeNode("f_int", [solve(expand2(g)), tree_form(var)])], bypart, sp, su)
    if eq is None:
        return None
    
    tmp2 = diffx(f, var) * eq
    eq2 = integratex(tmp2, depth-1, var, data+[TreeNode("f_int", [solve(expand2(tmp2)), tree_form(var)])], bypart, sp, su)
    
    if eq2 is None:
        return None
    
    return f * eq - eq2
tab = 0

def handle_intrec(special, equation):
    def extract_int(special, equation):
        
        if equation.name == "f_int" and equation == special:
            return parser.take_input("d")
        return TreeNode(equation.name, [extract_int(special, child) for child in equation.children])
    if equation is None:
        return None
    if "f_int" in str_form(equation):
        eq2 = extract_int(*copy.deepcopy([special, equation]))
        
        if "f_int" not in str_form(eq2):
            eq2 = solve(expand2(eq2))
            eq3 = solve(expand2(eq2-parser.take_input("d")))
            if eq3 != tree_form("d_0"):
                return inverse2(eq3, str_form(parser.take_input("d")))
        #return None
    return equation
def integratex(equation, depth, var, data, bypart, sp, su):
    if var not in ["v_0", "v_1", "v_2", "v_3"]:
        return None
    if depth <= 0:
        return None

    special = TreeNode("f_int", [solve(expand2(equation)), tree_form(var)])
    special2 = str_form(special)
    if sp and special2 in [str_form(item) for item in data[:-1]]:
        return special
    
    global tab
    plog(f"{'  '*tab}integrating {str(equation)} wrt {tree_form(var)}")
    
    if not contain(equation, tree_form(var)):
        return tree_form(var)*equation
    def varproc(term, var):
        if diffx(term, var).name[:2] == "d_":
            return diffx(term, var)
        return None
    
    equation = solve(quadratic(equation, True))
    
    equation = decompose(equation)
    plog(f"{'  '*tab}rewriting as {str(equation)}")
    
    eq2 = copy.deepcopy(equation)
    const = [x for x in factorgen(eq2) if var not in str_form(x)]
    tmp = tree_form("d_1")
    for item in const:
        tmp = tmp * item
    const = tmp
    const = solve(const)
    if const != tree_form("d_1"):
        tab += 1
        plog(f"{'  '*tab}taking the constant outside of the integral")
        tmp4 = solve(equation/const)
        r = integratex(tmp4, depth-1, var, data+[special], bypart, sp, su)
        
        plog(f"{'  '*tab}the solution is {str(r)}")
        tab -= 1
        if r is not None:
            return const*r
    orig = copy.deepcopy(equation)
    eq3 = solve(dowhile(replace_eqint(equation), expand2))
    eq2 = solve(dowhile(replace_eq3(equation), expand2))
    lst = [eq3, eq2]
    if solve(lst[1] - lst[0]) == tree_form("d_0"):
        lst = [lst[0]]
    for eq in lst:
        plog(f"{'  '*tab}rewriting as {str(eq)}")
        if eq.name == "f_add":
            
            success = True
            tab += 1
            plog(f"{'  '*tab}integating over sums")
            s = integratex(eq.children[0], depth - 1, var, data+[special], bypart, sp, su)
            s = handle_intrec(TreeNode("f_int", [solve(expand2(eq.children[0])), tree_form(var)]), s)
            plog(f"{'  '*tab}the solution is {str(s)}")
            tab -= 1
            if s is not None:
                for child in eq.children[1:]:
                    tab += 1
                    plog(f"{'  '*tab}integating over sums")
                    u = integratex(child, depth - 1, var, data+[special], bypart, sp, su)
                    u = handle_intrec(TreeNode("f_int", [solve(expand2(child)), tree_form(var)]), u)
                    plog(f"{'  '*tab}the solution is {str(u)}")
                    tab -= 1
                    if u is None:
                        
                        success = False
                        break

                    s += u
                if success:
                    return s
    eq = copy.deepcopy(orig)
    if eq.name == var:
        return tree_form(var)**tree_form("d_2")/tree_form("d_2")
    if eq.name == "f_pow" and eq.children[0].name == "s_e":
        tmp = varproc(eq.children[1], var)
        if tmp is not None:
            return eq/tmp
    if eq.name == "f_log" and eq.children[0].name == var:
        
        return eq*tree_form(var)-tree_form(var)
    tmp4 = structure(eq, tree_form('f_pow\n u_0\n u_1'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            expo = solve(tree_form(tmp4["u_1"])+tree_form("d_1"))
            base = tree_form(tmp4["u_0"])
            if expo == tree_form("d_0"):
                return base.fx("abs").fx("log")/tmp
            
            return base**expo / (expo*tmp)
    if eq.name == "f_cos":
        tmp = varproc(eq.children[0], var)
        if tmp is not None:
            
            return eq.children[0].fx("sin")/tmp
    if eq.name == "f_sin":
        tmp = varproc(eq.children[0], var)
        if tmp is not None:
            
            return tree_form("d_-1")*eq.children[0].fx("cos")/tmp
    tmp = structure(eq, tree_form('f_pow\n f_add\n  f_pow\n   u_0\n   d_2\n  p_0\n d_-1'))
    if tmp is not None:
        x = tree_form(tmp["u_0"])
        a = tree_form(tmp["p_0"])**(tree_form("d_2")**tree_form("d_-1"))
        if x.name == var:
            return (x/a).fx("arctan")/a
    if eq.name == "f_pow" and eq.children[0].name == "f_cos" and eq.children[1].name == "d_-2":
        tmp = varproc(eq.children[0].children[0], var)
        if tmp is not None:
            return eq.children[0].children[0].fx("tan")/tmp
    if eq.name == "f_pow" and eq.children[0].name == "f_sin" and eq.children[1].name == "d_-2":
        tmp = varproc(eq.children[0].children[0], var)
        if tmp is not None:
            return tree_form("d_-1")/(eq.children[0].children[0].fx("tan")*tmp)
    tmp4 = structure(eq, tree_form('f_mul\n f_sin\n  u_0\n f_pow\n  f_cos\n   u_0\n  d_-2'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form("d_1")/(tree_form(tmp4["u_0"]).fx("cos")*tmp)
    tmp4 = structure(eq, tree_form('f_mul\n f_cos\n  u_0\n f_pow\n  f_sin\n   u_0\n  d_-2'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            
            return tree_form("d_-1")/(tree_form(tmp4["u_0"]).fx("sin")*tmp)
    tmp4 = structure(eq, tree_form('f_mul\n f_pow\n  f_cos\n   u_0\n  d_-1\n f_sin\n  u_0'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form("d_-1")*tree_form(tmp4["u_0"]).fx("cos").fx("abs").fx("log")/tmp
    tmp4 = structure(eq, tree_form('f_mul\n f_pow\n  f_sin\n   u_0\n  d_-1\n f_cos\n  u_0'))
    if tmp4 is not None:
        tmp = varproc(tree_form(tmp4["u_0"]), var)
        if tmp is not None:
            return tree_form(tmp4["u_0"]).fx("sin").fx("abs").fx("log")/tmp
    tmp3 = structure(eq, tree_form('f_mul\n u_0\n f_pow\n  f_add\n   u_1\n   f_pow\n    u_2\n    u_3\n  f_mul\n   u_4\n   f_pow\n    u_5\n    d_-1'))
    if tmp3 is not None:
        tmp = [tree_form(tmp3[x]) for x in sorted(tmp3.keys(), key=lambda x: x)]
        if tmp[5] == tmp[3]:
            tmp6 = integratex(tmp[2]**tmp[4]*tmp[0]*(tmp[1]+tmp[2]**(tree_form("d_-1")*tmp[3]))**(tmp[4]/tmp[3]),depth - 1,  var, data+[special], bypart, sp, su)
            return tmp6

    if bypart:
        if eq.name == "f_mul" and len(eq.children) == 2:
            f, g = copy.deepcopy(eq.children)
            tab += 1
            plog(f"{'  '*tab}performing integration by parts f(x)={f} and g(x)={g}")
            tmp4 = byparts(f, g, depth, var, data+[special], bypart, sp, su)
            tab -= 1
            plog(f"{'  '*tab}the solution is {str(tmp4)}")
            if tmp4 is not None:
                
                return tmp4
            f, g = g, f
            tab += 1
            plog(f"{'  '*tab}performing integration by parts f(x)={f} and g(x)={g}")
            
            tmp4 = byparts(f, g, depth, var, data+[special], bypart, sp, su)
                
            tab -= 1
            plog(f"{'  '*tab}the solution is {str(tmp4)}")
            if tmp4 is not None:
                
                return tmp4
    
    if su:
        v2 = "v_"+str(int(var[2:])+1)
        for item in collec(eq, var):
            x = tree_form(v2)-item
            tab+= 1
            plog(f"{'  '*tab}substituting {str(item)}={tree_form(v2)}")
            tmp3 = subs1(equation, x, var, v2, depth, data+[special], bypart, sp, su)
            plog(f"{'  '*tab}the solution is {str(tmp3)}")
            tab -= 1
            if tmp3 is not None:
                return tmp3

    if not bypart:
        if eq.name == "f_mul" and len(eq.children) == 2:
            f, g = copy.deepcopy(eq.children)
            tab += 1
            plog(f"{'  '*tab}performing integration by parts f(x)={f} and g(x)={g}")
            tmp4 = byparts(f, g, depth, var, data+[special], bypart, sp, su)
            tab -= 1
            plog(f"{'  '*tab}the solution is {str(tmp4)}")
            if tmp4 is not None:
                
                return tmp4
            f, g = g, f
            tab += 1
            plog(f"{'  '*tab}performing integration by parts f(x)={f} and g(x)={g}")
            tmp4 = byparts(f, g, depth, var, data+[special], bypart, sp, su)
            tab -= 1
            plog(f"{'  '*tab}the solution is {str(tmp4)}")
            if tmp4 is not None:
                
                return tmp4
    return None
def formula_3(eq):
    if eq.name == "f_cos":
        return (parser.take_input("pi/2")-eq.children[0]).fx("sin")
    return TreeNode(eq.name, [formula_3(child) for child in eq.children])
def formula_4(eq):
    if eq.name == "f_add":
        for item in itertools.combinations(range(len(eq.children)), 2):
            if all(eq.children[item2].name == "f_sin" for item2 in item):
                a, b = eq.children[item[0]].children[0], eq.children[item[1]].children[0]
                rest = [item2 for index, item2 in enumerate(eq.children) if index not in item]
                if len(rest)==0:
                    rest = tree_form("d_0")
                else:
                    rest = summation(rest)
                two = tree_form("d_2")
                return rest + two*((a+b)/two).fx("sin")*((a-b)/two).fx("cos")
    return TreeNode(eq.name, [formula_4(child) for child in eq.children])
def formula_1(eq):
    if eq.name == "f_pow" and eq.children[1].name[:2] == "d_" and abs(int(eq.children[1].name[2:])) == 2 and eq.children[0].name == "f_cos":
        x = TreeNode("f_sin", [eq.children[0].children[0]])
        if int(eq.children[1].name[2:]) == 2:
            x = tree_form("d_1")-x*x
        else:
            x = tree_form("d_1")/(tree_form("d_1")-x*x)
        return x
    elif eq.name == "f_pow" and eq.children[1].name[:2] == "d_" and abs(int(eq.children[1].name[2:])) == 3 and eq.children[0].name == "f_cos":
        x = TreeNode("f_sin", [eq.children[0].children[0]])
        y = TreeNode("f_cos", [eq.children[0].children[0]])
        if int(eq.children[1].name[2:]) == 3:
            x = y-x*x*y
        else:
            x = tree_form("d_1")/(y-x*x*y)
        return x
    return TreeNode(eq.name, [formula_1(child) for child in eq.children])
def approx_limit(equation, var):
    orig = equation
    equation = approx(equation, var)
    if orig == equation:
        return None
    return equation
def find_limit(equation, depth=6, var=tree_form("v_0")):
    global tab
    ans = replace(equation, var, tree_form("d_0"))
    if depth == 0:
        return None
    try:
        compute(ans)
            
        plog(f"{'  '*tab}successful substitution")
        return expand4(solve(ans))
    except:
        pass
    tech = {0:"lhospital rule", 1:"converting into addition subproblems", 2: "approximation"}
    equation = solve(equation)
    
    equation2 = approx_limit(equation, var)
    if equation2 is not None:
    
        plog(f"{'  '*tab}after approximation limit {var}->0 {equation2}")
        tab += 1
        equation2 = find_limit(equation2, depth-1, var)
        tab -= 1
        if equation2 is not None:
            if len(varlist(equation2)) == 0 or var.name not in varlist(equation2):
                plog(f"{'  '*tab}solution is {equation2}")
                return equation2
        
    plog(f"{'  '*tab}trying to use {tech[0]}")
    ans = lhospital(equation, var)
    if ans is not None:
        
        plog(f"{'  '*tab}rewriting as limit {var}->0 {ans}")
        tab += 1
        equation2 = find_limit(copy.deepcopy(ans), depth-1, var)
        tab -= 1
        if equation2 is not None:
            if len(varlist(equation2)) == 0 or var.name not in varlist(equation2):
                plog(f"{'  '*tab}solution is {equation2}")
                return equation2

    equation2 = expand4(solve(equation))
    
    if equation2.name == "f_add":
        done = True
        s = tree_form("d_0")
        plog(f"{'  '*tab}applying limit across sums")
        for child in equation2.children:
            tab += 1
            out = find_limit(child, depth-1, var)
            tab -= 1
            if out is None:
                done = False
                break
            s+= out
        if done:
            s = expand4(solve(s))
            return s
    
    return None
def return_only_var(eq):
    output = []
    output2 = []
    def helper(eq, curse=False):
        nonlocal output
        if not curse and eq.name[:2] == "v_" and eq.name[2:].isdigit():
            output.append(eq.name)
        if eq.name[:2] == "v_" and eq.name[2:].isdigit():
            output2.append(eq.name)
        for child in eq.children:
            if child.name == "f_sin" or curse:
                helper(child, True)
            else:
                helper(child, False)
    helper(eq)
    output = list(set(output))
    output2 = list(set(output2))
    if len(output) == 0:
        output = [output2[0]]
    return output[0], None if len(output2)==1 else output[0]
def const_mul(equation):
    output = tree_form("d_1")
    for item in factorgen(equation):
        if "v_" not in str_form(item):
            output = output * item
    return output
def degree(eq, var):
    if any("f_" in x in str_form(x) for x in "sin cos tan".split(" ")):
        return 4
    count = 0
    while True:
        if eq == tree_form("d_0") or count==4:
            break
        eq = diffx(eq, var)
        count += 1
    return count
def equal(a, b):
    x = a-b
    x = expand2(x)
    x = solve(x)
    if x == tree_form("d_0"):
        return True
    return False
def replace2(equation, find, r):
    if equation.name == "f_int":
        return equation
    if equal(equation, find):
        return r
    col = TreeNode(equation.name, [])
    for child in equation.children:
        col.children.append(replace2(child, find, r))
    return col
def contain(eq, term):
    if equal(eq, term):
        return True
    for child in eq.children:
        if contain(child, term):
            return True
    return False
def clearv(eq, varname=None):
    def vlist(eq):
        out = []
        if eq.name[:2] == "v_":
            out.append(eq.name)
        for child in eq.children:
            out += vlist(child)
        return sorted(list(set(out)), key=lambda x: int(x[2:]))
    if eq.name == "f_eq":
        equation=eq.children[0]
        arr = None
        if varname is None:
            arr = [x for x in factorgen(equation) if "v_" in str_form(x)]
        else:
            arr = [x for x in factorgen(equation) if varname in vlist(x)]
        for i in range(len(arr)-1,-1,-1):
            if arr[i].name == "f_pow" and arr[i].children[1].name[:2] == "d_" and int(arr[i].children[1].name[2:]) < 0:
                arr.pop(i)
        m = arr[0]
        for item in arr[1:]:
            m = m * item
        return TreeNode(eq.name, [m, tree_form("d_0")])
    return TreeNode(eq.name, [clearv(child) for child in eq.children])
def fraction2(equation):
    #equation = formula_1(expand(common2(equation)))
    equation = expand(common2(equation))
    equation = solve(equation)
    return equation
def inverse(rhs,term):
    lhs = tree_form("d_0")
    count = 6
    while not equal(rhs, term):
        if rhs.name == "f_add":
            if all(term in factorgen(child) for child in rhs.children):
                newrhs = solve(expand2(rhs/term))
                if not contain(newrhs, term):
                    rhs = term * newrhs
            else:
                for i in range(len(rhs.children)-1,-1,-1):
                    if not contain(rhs.children[i], term):
                        lhs = lhs - rhs.children[i]
                        rhs.children.pop(i)
        elif rhs.name == "f_mul":
            for i in range(len(rhs.children)-1,-1,-1):
                if not contain(rhs.children[i], term):
                    lhs = lhs / rhs.children[i]
                    rhs.children.pop(i)
        elif rhs.name == "f_pow" and contain(rhs.children[0], term):
            lhs = lhs ** (tree_form("d_1")/rhs.children[1])
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_sin" and contain(rhs.children[0], term):
            lhs = lhs.fx("arcsin")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_arcsin" and contain(rhs.children[0], term):
            lhs = lhs.fx("sin")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_arccos" and contain(rhs.children[0], term):
            lhs = lhs.fx("cos")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_cos" and contain(rhs.children[0], term):
            lhs = lhs.fx("arccos")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_log" and contain(rhs.children[0], term):
            lhs = tree_form("s_e")**lhs
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_pow" and rhs.children[0].name == "s_e" and contain(rhs.children[1], term):
            lhs = lhs.fx("log")
            rhs = copy.deepcopy(rhs.children[1].fx("log"))
        elif rhs.name == "f_tan" and contain(rhs.children[0], term):
            lhs = lhs.fx("arctan")
            rhs = copy.deepcopy(rhs.children[0])
        elif rhs.name == "f_arctan" and contain(rhs.children[0], term):
            lhs = lhs.fx("tan")
            rhs = copy.deepcopy(rhs.children[0])
        if len(rhs.children) == 1:
            rhs = rhs.children[0]
        count -= 1
        if count == 0:
            return None
    return solve(lhs)
def inverse2(eq, var="v_0"):
    def newfx(eq):
        eq = expand(common2(eq))
        eq = clearv(TreeNode("f_eq", [eq, tree_form("d_0")])).children[0]
        eq = solve(eq)
        return eq
    eq = dowhile(eq, newfx)
    eq2 = inverse(eq, tree_form(var))
    if eq2 is None:
        return None
    tmp = solve(expand2(eq2))
    return tmp
def subs1(equation, term, v1, v2, depth, data, bypart, sp, su):

    origv2 = copy.deepcopy(v2)
    equation = solve(equation)
    eq = equation
    termeq = term
    t = inverse2(copy.deepcopy(termeq), v1)
    g = inverse2(termeq, v2)
    
    if g is None:
        return None
    if t is None:
        return None
        eq2 = diffx(g, v1)
        equation = solve(eq/eq2)
        equation = replace2(equation, g, tree_form(v2))
        if v1 in str_form(equation):
            return None
    else:
        t = expand2(t)
        eq = replace2(eq, tree_form(v1), t)
        
        
        eq2 = replace2(diffx(g, v1), tree_form(v1), t)
        equation = eq/eq2
        equation = solve(equation)
        
    lst = [ equation]
    for eq in lst:
        if v1 in str_form(eq):
            continue
        
        eq = solve(dowhile(eq, expand2))
        
        tmp = integratex(eq, depth - 1, origv2, data, bypart, sp, su)
        if tmp is None:
            continue
        tmp = replace2(tmp, tree_form(v2), g)
        return tmp
    return None
def rref(matrix):
    rows, cols = len(matrix), len(matrix[0])
    lead = 0
    for r in range(rows):
        if lead >= cols:
            return matrix
        i = r
        while solve(matrix[i][lead]) == tree_form("d_0"):
            i += 1
            if i == rows:
                i = r
                lead += 1
                if lead == cols:
                    return matrix
        matrix[i], matrix[r] = matrix[r], matrix[i]
        lv = matrix[r][lead]
        matrix[r] = [solve(m / lv) for m in matrix[r]]
        for i in range(rows):
            if i != r:
                lv = matrix[i][lead]
                matrix[i] = [solve(m - lv * n) for m, n in zip(matrix[i], matrix[r])]
        lead += 1
    return matrix
def islinear(eq, fxconst):
    eq =solve(eq)
    if eq.name == "f_pow" and fxconst(eq):#"v_" in str_form(eq):
        return False
    for child in eq.children:
        out = islinear(child, fxconst)
        if not out:
            return out
    return True
def linear(eqlist, fxconst):
    final = []
    extra = []
    for i in range(len(eqlist)-1,-1,-1):
        if eqlist[i].name == "f_mul" and not islinear(expand2(eqlist[i]), fxconst):
            if "v_" in str_form(eqlist[i]):
                eqlist[i] = TreeNode("f_mul", [child for child in eqlist[i].children if fxconst(child)])
            if all(islinear(child, fxconst) for child in eqlist[i].children):
                for child in eqlist[i].children:
                    extra.append(TreeNode("f_eq", [child, tree_form("d_0")]))
                eqlist.pop(i)
            else:
                final.append(TreeNode("f_eq", [eqlist[i], tree_form("d_0")]))
                eqlist.pop(i)
    
    if extra != []:
        final.append(TreeNode("f_or", extra))
    if eqlist == []:
        if len(final)==1:
            
            return final[0]
        return TreeNode("f_and", final)
    eqlist = [eq for eq in eqlist if fxconst(eq)]
    if not all(islinear(eq, fxconst) for eq in eqlist):
        return TreeNode("f_and", copy.deepcopy(final+eqlist))
    vl = []
    def varlist(eq, fxconst):
        nonlocal vl
        if eq.name[:2] == "v_" and fxconst(eq):
            vl.append(eq.name)
        for child in eq.children:
            varlist(child, fxconst)
    for eq in eqlist:
        varlist(eq, fxconst)
    vl = list(set(vl))
    if len(vl) > len(eqlist):
        return TreeNode("f_and", final+[TreeNode("f_eq", [x, tree_form("d_0")]) for x in eqlist])
    m = []
    for eq in eqlist:
        s = copy.deepcopy(eq)
        row = []
        for v in vl:
            row.append(diffx(eq, v))
            s = replace(s, tree_form(v), tree_form("d_0"))
        row.append(s)
        m.append(row)
    for i in range(len(m)):
        for j in range(len(m[i])):
            m[i][j] = solve(expand2(m[i][j]))
     
    m = rref(m)
    for i in range(len(m)):
        for j in range(len(m[i])):
            m[i][j] = fraction2(fraction2(m[i][j]))

    for item in m:
        if all(item2==tree_form("d_0") for item2 in item[:-1]) and item[-1] != tree_form("d_0"):
            return tree_form("d_false")
    
    output = []
    for index, row in enumerate(m):
        count = 0
        for item in row[:-1]:
            if item == tree_form("d_1"):
                count += 1
                if count == 2:
                    break
            elif item == tree_form("d_0") and count == 1:
                break
        if count == 0:
            continue
        output.append(tree_form(vl[index])+row[-1])
    if len(output) == 1 and len(final)==0:
        return TreeNode("f_eq", [output[0], tree_form("d_0")])
    return TreeNode("f_and", final+[TreeNode("f_eq", [x, tree_form("d_0")]) for x in output])


def quadratic(equation, integration=False, to_compute=None, root=False, com = False):
    output = None
    def quad(eq):
        nonlocal to_compute
        nonlocal root
        if to_compute is None:
            tmp, to_compute = return_only_var(eq)
        else:
            tmp = tree_form(to_compute)
        if tmp is not None:
            eq = expand2(eq)
            eq = solve(eq)
            
            if eq.name == "f_add":
                
                for i in range(len(eq.children)):
                    if "v_" in str_form(eq.children[i]) and eq.children[i].name != "f_mul":
                        eq.children[i] = tree_form("d_1")*eq.children[i]
                if all("v_" not in str_form(child) for child in eq.children):
                    eq = eq + tree_form("d_0")
                tmp = structure(eq, tree_form('f_add\n f_mul\n  f_pow\n   u_0\n   d_2\n  p_0\n f_mul\n  u_0\n  p_1\n p_2'), {}, to_compute)
                if tmp is not None:
                    
                    a, b, c = [tree_form(x) for x in [tmp["p_0"], tmp["p_1"], tmp["p_2"]]]
                    var = tree_form(tmp["u_0"])
                    d = (b**tree_form("d_2") - tree_form("d_4")*a*c)
                    d = solve(expand2(d))
                    if compute(d)>=0 or to_compute is not None or com:
                        r1 = ( tree_form("d_-1")*b + (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r2 = ( tree_form("d_-1")*b - (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r1, r2 = solve(r1), solve(r2)
                        if root:
                            
                            return [inverse2(var-r1, varlist(var)[0]), inverse2(var-r2, varlist(var)[0])]
                        return a*solve((var-r1)*(var-r2))
                    elif integration:
                        r1 = b/(tree_form("d_2")*a)
                        r2 = c/a - b**tree_form("d_2")/(tree_form("d_4")*a**tree_form("d_2"))
                        return a*solve((var+r1)**tree_form("d_2") + r2)
                
                tmp = structure(eq, tree_form('f_add\n f_mul\n  p_0\n  f_pow\n   u_0\n   d_2\n f_mul\n  p_1\n  u_0'), {}, to_compute)
                if tmp is not None:
                    a, b, c = [tree_form(x) for x in [tmp["p_0"], tmp["p_1"], "d_0"]]
                    var = tree_form(tmp["u_0"])
                    d = (b**tree_form("d_2") - tree_form("d_4")*a*c)
                    d = solve(expand2(d))
                    if compute(d)>=0 or to_compute is not None or com:
                        r1 = ( tree_form("d_-1")*b + (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r2 = ( tree_form("d_-1")*b - (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r1, r2 = solve(r1), solve(r2)
                        if root:
                            return [inverse2(var-r1, varlist(var)[0]), inverse2(var-r2, varlist(var)[0])]
                        return a*solve((var-r1)*(var-r2))
                    elif integration:
                        r1 = b/(tree_form("d_2")*a)
                        r2 = c/a - b**tree_form("d_2")/(tree_form("d_4")*a**tree_form("d_2"))
                        if root:
                            
                            return [inverse2(var-r1, varlist(var)[0]), inverse2(var-r2, varlist(var)[0])]
                        return a*solve((var+r1)**tree_form("d_2") + r2)
                    return a*solve(var**tree_form("d_2") + c/a)
                
                tmp = structure(eq, tree_form('f_add\n f_mul\n  f_pow\n   u_0\n   d_2\n  p_0\n p_1'), {}, to_compute)
                if tmp is not None:
                    a, b, c = [tree_form(x) for x in [tmp["p_0"], "d_0", tmp["p_1"]]]
                    var = tree_form(tmp["u_0"])
                    d = (b**tree_form("d_2") - tree_form("d_4")*a*c)
                    d = solve(expand2(d))
                    if compute(d)>=0 or to_compute is not None or com:
                        r1 = ( tree_form("d_-1")*b + (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r2 = ( tree_form("d_-1")*b - (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r1, r2 = solve(r1), solve(r2)
                        if root:
                            return [inverse2(var-r1, varlist(var)[0]), inverse2(var-r2, varlist(var)[0])]
                        return a*solve((var-r1)*(var-r2))
                    elif integration:
                        r1 = b/(tree_form("d_2")*a)
                        r2 = c/a - b**tree_form("d_2")/(tree_form("d_4")*a**tree_form("d_2"))
                        return a*solve((var+r1)**tree_form("d_2") + r2)
                    return a*solve(var**tree_form("d_2") + c/a)
                tmp = structure(eq, tree_form('f_add\n f_mul\n  u_0\n  p_0\n p_1'), {}, to_compute)
                if tmp is not None:
                    
                    a, b, c = [tree_form(x) for x in ["d_0", tmp["p_0"], tmp["p_1"]]]
                    var = tree_form(tmp["u_0"])
                    r1= b*solve(var+c/b)
                    if root:
                        return [r1]
                    return r1
                tmp = structure(eq, tree_form('f_add\n f_add\n  f_mul\n   p_0\n   f_pow\n    u_0\n    d_4\n  f_mul\n   p_1\n   f_pow\n    u_0\n    d_2\n p_2'), {}, to_compute)
                if tmp is not None:
                    
                    a, b, c = [tree_form(x) for x in [tmp["p_0"], tmp["p_1"], tmp["p_2"]]]
                    var = tree_form(tmp["u_0"])
                    d = (b**tree_form("d_2") - tree_form("d_4")*a*c)
                    d = solve(expand2(d))
                    if compute(d)>=0 or to_compute is not None or com:
                        r1 = ( tree_form("d_-1")*b + (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        r2 = ( tree_form("d_-1")*b - (b**tree_form("d_2") - tree_form("d_4")*a*c)** (tree_form("d_1")/tree_form("d_2") ) )/(tree_form("d_2")*a)
                        
                        r1, r2 = r1**parser.take_input("1/2"), r2**parser.take_input("1/2")
                        r1, r2 = solve(r1), solve(r2)
                        r3, r4 = -r1, -r2
                        
                        if root:
                            
                            return [inverse2(var-r1, varlist(var)[0]), inverse2(var-r2, varlist(var)[0]), inverse2(var-r3, varlist(var)[0]), inverse2(var-r4, varlist(var)[0])]
                        return a*solve((var-r1)*(var-r2))
        return None
    if equation.name == "f_add":
        tmp = quad(copy.deepcopy(equation))
        if tmp is not None:
            if not isinstance(tmp, TreeNode):
                if root:
                    return tmp
                tmp = tmp[0] * tmp[1]
            return tmp
    if equation.name == "f_eq" and root:
        return quadratic(equation.children[0], integration, to_compute, root)
    return TreeNode(equation.name, [quadratic(child, integration, to_compute, root) for child in equation.children])
def zero(eq):
    eq = copy.deepcopy(eq)
    if eq.name != "f_list":
        return tree_form("d_0")
    return TreeNode(eq.name, [zero(child) for child in eq.children])
def eqhandle(eq):
    if eq.name in ["f_eq", "f_lt", "f_gt"]:
        if "f_list" in str_form(eq):
            if eq.children[1] != zero(eq.children[1]):
                eq.children[0] = eqhandle(eq.children[1] - eq.children[0])
                eq.children[1] = zero(eq.children[1])
            else:
                eq.children[0] = eqhandle(eq.children[0])
            return eq
        return TreeNode(eq.name, [eqhandle(eq.children[1]-eq.children[0]), zero(eq.children[1])])
    return TreeNode(eq.name, [eqhandle(copy.deepcopy(child)) for child in eq.children])
def rmeq(eq):
    if eq.name == "f_eq":
        return rmeq(eq.children[0])
    return TreeNode(eq.name, [rmeq(child) for child in eq.children])
def find_abs(eq, eqlst):
    out = []
    if eq.name == "f_abs":
        for item in eqlst:
            if solve(expand2(item - eq.children[0])) == tree_form("d_0"):
                return solve(expand2(eq.children[0] * tree_form("d_-1")))
            elif solve(expand2(item + eq.children[0])) == tree_form("d_0"):
                
                return solve(expand2(eq.children[0]))
    return TreeNode(eq.name, [find_abs(child, eqlst) for child in eq.children])

def mat0(eq, lst=None):
    def findeq(eq):
        out = []
        if "f_list" not in str_form(eq) and "f_eq" not in str_form(eq):
            return [str_form(eq)]
        else:
            for child in eq.children:
                out += findeq(child)
        return out
    eqlist = findeq(eq)
    eqlist = [tree_form(x) for x in eqlist]
    eqlist = [rmeq(x) for x in eqlist]
    eqlist = [TreeNode("f_mul", factorgen(x)) for x in eqlist if x != tree_form("d_0")]
    eqlist = [x.children[0] if len(x.children) == 1 else x for x in eqlist]
    out = None
    
    if lst is None:
        out = linear(copy.deepcopy(eqlist), lambda x: "v_" in str_form(x))
    else:
        out = linear(copy.deepcopy(eqlist), lambda x: any(contain(x, item) for item in lst))
    def rms(eq):
        if eq.name in ["f_and", "f_or"] and len(eq.children) == 1:
            return eq.children[0]
        return TreeNode(eq.name, [rms(child) for child in eq.children])
    return rms(out)
def and0(eq, lst=None):
    if eq.name == "f_and":
        eq2 = copy.deepcopy(eq)
        eq2.name = "f_list"
        return mat0(eq2, lst)
    elif eq.name == "f_eq":
        return mat0(eq, lst)
    return TreeNode(eq.name, [and0(child, lst) for child in eq.children])


def mat00(eq):
    if eq.name == "f_eq":
        if "f_list" in str_form(eq):
            return mat0(eq.children[0])
        return TreeNode(eq.name, [eq.children[0], tree_form("d_0")])
    return TreeNode(eq.name, [mat00(child) for child in eq.children])

def extract(eq, s, lst):
    lst2 = [parser.take_input(x).name for x in lst]
    out = {}
    out2 = TreeNode("f_and", [])
    if eq.name == "f_and":
        for child in eq.children:
            v = list(set(lst2) & set(varlist(child)))
            if v == []:
                out2.children.append(child)
                continue
            else:
                v = v[0]
            n = solve(child.children[0] - tree_form(v))
            out[v] = n
    def r(s):
        nonlocal out
        if s.name in out.keys():
            return out[s.name]
        return TreeNode(s.name, [r(child) for child in s.children])
    return r(s), out2
def normal(string, variable=None):
    equation = None
    if variable is not None:
        equation = parser.take_input(string,list(variable.keys()))
    else:
        equation = parser.take_input(string,None)
    def req(eq, key, variable, vl):
        if eq.name == "f_"+key:
            tmp = variable[key]
            for i in range(len(eq.children)):
                tmp = replace(copy.deepcopy(tmp), tree_form(vl[i]), eq.children[i])
            return tmp
        return TreeNode(eq.name, [req(child, key, variable, vl) for child in eq.children])
    if variable is not None:
        for key in variable.keys():
            equation = req(equation, key, variable, varlist(variable[key]))
    alloclst = []
    for i in range(0,26):
        if "v_"+str(i) not in varlist(equation):
            alloclst.append("v_"+str(i))
    va, vb, vc, vd = [str(tree_form(item)) for item in alloclst[:4]]
    linesegment = parser.take_input(f"[[{va},{vb}],[{vc},{vd}]]")
    edit = False
    if "f_linesegment" in str_form(equation):
        edit = True
        eq3 = parser.take_input("linesegment(0)")
        equation = replace(equation, eq3, linesegment)
    equation = matrix_simp2(equation)
    equation = eqhandle(equation)
    if "f_list" in str_form(equation):
        equation = copy.deepcopy(convert_sub2neg(equation))
        equation = matrix_simp(equation)
        equation = matrix_simp(equation)
    if edit:
        eq2, eq3 = extract(and0(equation), linesegment, [va, vb, vc, vd])
        equation = TreeNode("f_eq", [parser.take_input("linesegment(0)"), eq2])
        if eq3.children:
            equation = TreeNode("f_and", eq3.children + [equation])
    equation = rmdeg(equation)
    return equation
def formula_2(eq, var="v_0"):
    orig = copy.deepcopy(eq)
    eq = replace(eq, tree_form(var).fx("sin")**tree_form("d_2"), (tree_form("d_1") - tree_form(var).fx("cos")**tree_form("d_2")))
    eq = flatten_tree(eq)
    eq = solve(expand2(eq))
    if "f_sin" not in str_form(eq):
        return eq
    eq = orig
    orig = copy.deepcopy(eq)
    eq = replace(eq, tree_form(var).fx("cos")**tree_form("d_2"), (tree_form("d_1") - tree_form(var).fx("sin")**tree_form("d_2")))
    eq = flatten_tree(eq)
    eq = solve(expand2(eq))
    if "f_cos" not in str_form(eq):
        return eq
    return eq
intconst = ["v_c"+str(i) for i in range(1,10)]
def allocvar():
    global intconst
    return tree_form(intconst.pop(0))
def intreg(eq):
    if eq.name == "f_integrate":
        v = allocvar()
        return integratex(eq.children[0], eq.children[1].name) + v
    return TreeNode(eq.name, [intreg(child) for child in eq.children])
def summation(lst):
    s = lst[0]
    for item in lst[1:]:
        s += item
    return s
def product(lst):
    if lst == []:
        return tree_form("d_1")
    s = lst[0]
    for item in lst[1:]:
        s *= item
    return s
def epowersplit(eq):
    if eq.name == "f_pow" and eq.children[1].name == "f_add":
        return product([eq.children[0]**child for child in eq.children[1].children])
    return TreeNode(eq.name, [epowersplit(child) for child in eq.children])
def esolve(s):
    if s.name == "f_add" and "f_log" in str_form(s):
        return product([tree_form("s_e")**child for child in s.children]) - tree_form("d_1")
    return TreeNode(s.name, [esolve(child) for child in s.children])
def diffsolve_sep2(eq):
    global tab
    s = []
    eq = solve4(expand2(eq))
    eq = e1(eq)
    if eq.name == "f_add" and all(len(varlist(child))==1 and [str_form(x) for x in factorgen(copy.deepcopy(child))].count(str_form(tree_form(varlist(child)[0]).fx("dif")))==1 for child in eq.children):
        for child in eq.children:
            v = varlist(child)[0]
            v2 = tree_form(v).fx("dif")
            child = replace(child, v2, tree_form("d_1"))
            child = solve(child)
            tab += 1
            
            tmp6 = copy.deepcopy(integratex(child, 4, v))
            s.append(tmp6)
            plog(f"{'  '*tab}the solution is {tmp6}")
            tab -= 1
            if s[-1] is None:
                return None
        s.append(allocvar())
    else:
        return None
    s = summation(s)
    s = solve4(e0(s))
    return groupe(s)
def e0(eq):
    return TreeNode("f_eq", [eq, tree_form("d_0")])
def e1(eq):
    if eq.name == "f_eq":
        eq = eq.children[0]
    return eq
def groupe(eq):
    eq = esolve(eq)
    eq = solve4(eq)
    eq = fraction2(eq)
    eq = clearv(eq)
    eq = epowersplit(eq)
    return eq
from collections import Counter
def multiset_intersection(*lists):
    counters = list(map(Counter, lists))
    common = counters[0]
    for c in counters[1:]:
        common = common & c
    return list(common.elements())
def subtract_sublist(full_list, sublist):
    c_full = Counter(full_list)
    c_sub = Counter(sublist)
    result = c_full - c_sub
    tmp = list(result.elements())
    if tmp == []:
        return [tree_form("d_1")]
    return tmp

def term_common2(eq):
    if eq.name != "f_add":
        return eq
    s = []
    arr = [factorgen(child) for child in eq.children]
    s = multiset_intersection(*arr)
    return product(s)*summation([product(subtract_sublist(factorgen(child), s)) for child in eq.children])
def term_common(eq):
    if eq.name == "f_add":
        return solve(term_common2(eq))
    return solve(product([term_common2(item) for item in factorgen(eq)]))

def take_common(eq):
    if eq.name == "f_add":
        eq = term_common(eq)
        if eq.name == "f_add":
            for i in range(len(eq.children)-1,1,-1):
                for item in itertools.combinations(range(len(eq.children)), i):
                    eq2 = summation([item2 for index, item2 in enumerate(eq.children) if index in item])
                    eq2 = term_common(eq2)
                    if eq2.name == "f_mul":
                        return take_common(solve(summation([item2 for index, item2 in enumerate(eq.children) if index not in item]) + eq2))
        return eq
    return term_common(eq)
def diffsolve_sep(eq):
    eq = groupe(eq)
    def inversediff(lhs, rhs):
        count = 4
        while contain(rhs, tree_form("v_1")) or contain(lhs, tree_form("v_0")):
            success = False
            if rhs.name == "f_add":
                for i in range(len(rhs.children)-1,-1,-1):
                    if not contain(rhs.children[i], tree_form("v_0")) or str_form(tree_form("v_1").fx("dif")) in [str_form(x) for x in factorgen(rhs.children[i])]:
                        if "v_" in str_form(rhs.children[i]):
                            success = True
                        lhs = lhs - rhs.children[i]
                        rhs.children.pop(i)
            elif rhs.name == "f_mul":
                for i in range(len(rhs.children)-1,-1,-1):
                    if not contain(rhs.children[i], tree_form("v_0")):
                        if "v_" in str_form(rhs.children[i]):
                            success = True
                        lhs = lhs / rhs.children[i]
                        rhs.children.pop(i)
            if len(rhs.children) == 1:
                rhs = rhs.children[0]
            rhs, lhs = copy.deepcopy([solve4(lhs), solve4(rhs)])
            if rhs.name == "f_add":
                for i in range(len(rhs.children)-1,-1,-1):
                    if not contain(rhs.children[i], tree_form("v_1")) or str_form(tree_form("v_0").fx("dif")) in [str_form(x) for x in factorgen(rhs.children[i])]:
                        if "v_" in str_form(rhs.children[i]):
                            success = True
                        lhs = lhs - rhs.children[i]
                        rhs.children.pop(i)
            elif rhs.name == "f_mul":
                for i in range(len(rhs.children)-1,-1,-1):
                    if not contain(rhs.children[i], tree_form("v_1")):
                        if "v_" in str_form(rhs.children[i]):
                            success = True
                        lhs = lhs / rhs.children[i]
                        rhs.children.pop(i)
            rhs, lhs = copy.deepcopy([solve4(lhs), solve4(rhs)])
            if not success:
                lhs, rhs = term_common(lhs),term_common(rhs)
            count -= 1
            if count == 0:
                return None
        return solve4(e0(lhs-rhs))
    plog(f"{'  '*tab}rewritting as {eq}")
    eq = inversediff(tree_form("d_0"), copy.deepcopy(eq.children[0]))
    plog(f"{'  '*tab}rearranged as {eq}")
    if eq is None:
        return None
    tmp = diffsolve_sep2(eq)
    if tmp is not None:
        return tmp
    return None
def diffsolve(eq):
    global tab
    plog(f"{'  '*tab}solving the differential equation {eq}")
    eq = epowersplit(eq)
    eq = fraction2(eq)
    
    eq = diffsolve_sep(eq)
    
    plog(f"{'  '*tab}the solution is {eq}")
    return eq

def lineardiff(eq):
    eq = epowersplit(eq)
    eq = fraction2(eq)
    eq = expand2(eq)
    eq = clearv(eq)
    eq = fraction2(eq)
    eq = clearv(eq)
    eq = eq.children[0]
    eq = TreeNode(eq.name, [term_common(eq), tree_form("d_0")])
    eq = solve(eq)
    eq = clearv(eq)
    eq = eq.children[0]
    if eq.name == "f_add" and len(eq.children) == 3:
        for item in itertools.permutations(eq.children):
            if parser.take_input("dif(y)") in factorgen(item[0]) and\
               parser.take_input("y") in factorgen(item[1]) and parser.take_input("dif(x)") in factorgen(item[1]) and\
               parser.take_input("dif(x)") in factorgen(item[2]):
                m = solve(item[0]/parser.take_input("dif(y)"))
                item2 = []
                item2.append(item[0])
                item2.append(solve(item[1]/m))
                item2.append(solve(item[2]/m))
                print(item2)
                p = solve(item2[1]/parser.take_input("y*dif(x)"))
                q = solve(item2[2]/parser.take_input("-dif(x)"))
                f = tree_form("s_e") ** integratex(p, 3)
                return TreeNode("f_eq", [parser.take_input("y")*f - integratex(q*f, 3) + allocvar(), tree_form("d_0")])
def varlist(eq):
    out = []
    if eq.name[:2] == "v_":
        out.append(eq.name)
    for child in eq.children:
        out += varlist(child)
    return sorted(list(set(out)), key=lambda x: int(x[2:]))
def subshomodif(eq):
    global tab
    def dodiff(eq):
        if eq.name == "f_dif" and eq.children[0].name[:2] != "v_":
            return diff(eq.children[0])
        return TreeNode(eq.name, [dodiff(child) for child in eq.children])
    eq1 = replace(eq, parser.take_input("y"), parser.take_input("x*z"))
    eq1 = solve(eq1)
    eq2 = replace(eq, parser.take_input("x"), parser.take_input("y*z"))
    eq2 = solve(eq2)
    plog(f"{'  '*tab}solving the differential equation {eq}")
    for index, item in enumerate([eq1, eq2]):
        tmp = varlist(item)
        if len(tmp) <=2:
             item = dodiff(item)
             if index == 0:
                 item = replace(item, parser.take_input("z"), parser.take_input("y"))
                 plog(f"{'  '*tab}homogenous, substituting z=y/x")
             else:
                 item = replace(item, parser.take_input("z"), parser.take_input("x"))
                 plog(f"{'  '*tab}homogenous, substituting z=x/y")
             item = expand2(item)
             item = solve(item)
             tab += 1
             item2 = diffsolve(item)
             
             tab -= 1
             if item2 is None:
                 continue
             if index == 0:
                 item2 = replace(item2, parser.take_input("y"), parser.take_input("z"))
                 item2 = replace(item2, parser.take_input("z"), parser.take_input("y/x"))
             else:
                 item2 = replace(item2, parser.take_input("x"), parser.take_input("z"))
                 item2 = replace(item2, parser.take_input("z"), parser.take_input("x/y"))
             return item2
    plog(f"{'  '*tab}not homogenous")
def homo(eq):
    eq = inverse2(eq.children[0], str_form(tree_form("v_1").fx("dif")))
    eq = solve(eq/parser.take_input("dif(x)"))
def expand3(equation):
    if equation.name == "f_eq":
        eq2 = copy.deepcopy(equation.children[1])
        equation = TreeNode("f_eq", [expand2(e1(equation)), eq2])
    else:
        equation = expand2(equation)
    equation = solve(equation)
    return equation
def expand4(equation):
    return solve(expand3(replace_eq(replace_eq2(equation))))
def is_valid(eq):
    try:
        val =compute(eq)
    except:
        return False
    return True
def trigcommand(equation):
    equation = solve(expand2(replace_eq(equation)))
    print(equation)
    return equation
def trigsolve(eq, compute_eq=None):
    lhs = eq
    rhs = parser.take_input("0")
    lhs = replace_eq2(lhs)
    lhs = dowhile(lhs, fraction2)
    lhs = clearv(lhs)
    lhs = expand3(lhs)
    lhs2 = replace(lhs, parser.take_input("sin(x)"), parser.take_input("y"))
    lhs2 = inverse2(lhs2.children[0], "v_1")
    lhs3 = replace(lhs2, parser.take_input("y"), parser.take_input("sin(x)"))

    lhs = lhs3**parser.take_input("2")-parser.take_input("1-cos(x)^2")
    
    lhs = expand3(solve(lhs))
    print(f"square of the rewritten expression = {str(lhs)}")
    lhs = list(quadratic(lhs, False, "v_0", True, False))
    
    if compute_eq is None:
        return lhs[0]
    lhs += [-x for x in lhs]
    root = []
    for item in lhs[::-1]:
        print(f"x can be {str(item)}")
        if not is_valid(item):
            continue
        ques = replace(lhs3  - parser.take_input("sin(x)"), parser.take_input("x"), item)
        ques = expand3(solve(ques))
        ques = trigcommand(ques)
        print(f"will be checking if {str(ques)} is 0")
        if ques != parser.take_input("0"):
            continue
        rhs = replace(compute_eq, parser.take_input("x"), item)
        rhs = replace_eq2(rhs)
        rhs = solve(rhs)
        if rhs is not None:
            rhs = expand3(solve(rhs))
            return rhs
    return None
    if len(root) == 1:
        return root[0]
    return TreeNode("f_and", root)
def contain_fx(eq, mode=False):
    def which(eq):
        output = []
        if eq.name[:2] == "f_" and len(eq.name) == 3:
            if mode and eq.children[0].name[:2] == "d_":
                output.append(str_form(eq))
            elif not mode and eq.children[0].name[:2] =="v_":
                output.append(str_form(eq))
            
        for child in eq.children:
            output += which(child)
        return list(set(output))
    output = []
    if eq.name == "f_eq":
        for w in which(eq.children[0]):
            w = tree_form(w)
            if not mode:
                out = inverse(eq.children[0], w)
                output.append([w.name, replace(out, w.children[0], parser.take_input("x"))])
            else:
                output.append(w)
    for child in eq.children:
        output += contain_fx(child)
    return output

def isfx(eq):
    if eq.name[:2] == "f_" and len(eq.name) == 3:
        return True
    if any(isfx(child) for child in eq.children):
        return True
    return False
def fxmaker(eqlst):
    fx_dic = {}

    for i in range(len(eqlst)):
        out = contain_fx(copy.deepcopy(eqlst[i]))
        if len(out) != 0:
            for item in out:
                fx_dic[item[0]] = item[1]
             
    for i in range(len(eqlst)):
        out = contain_fx(eqlst[i], True)
        if len(out) != 0:
            for item in out:
                
                ans = replace(fx_dic[item.name], parser.take_input("x"), item.children[0])
                
                eqlst[i] = copy.deepcopy(replace(eqlst[i], item, ans))
    return eqlst, fx_dic
def cs(eqlst, haveto):
    eqlst, fx_dic = fxmaker(eqlst)
    def fxfx(eq):
        nonlocal fx_dic
        if eq.name in fx_dic.keys():
            return replace(fx_dic[eq.name], parser.take_input("x"), eq.children[0])
        return TreeNode(eq.name, [fxfx(child) for child in eq.children])
    return expand3(solve(fxfx(haveto)))
def find(eqlst, to_compute):
    
    eqlst, fx_dic = fxmaker(eqlst)
    for i in range(len(eqlst)-1, -1, -1):
        if isfx(eqlst[i]):
            
            eqlst.pop(i)
            continue
        v = list(set(varlist(eqlst[i])) & set([str_form(item) for item in to_compute]))
        if len(v) != 1:
            continue
        eqlst[i] = copy.deepcopy(expand3(solve(eqlst[i])))
        out = quadratic(eqlst[i], False, v[0], True, False)
        if isinstance(out, TreeNode):
            continue
        new = list(out)
        if isinstance(new, list):
            new= TreeNode("f_or", [TreeNode("f_eq", [tree_form(v[0])-item, tree_form("d_0")]) for item in new])
            eqlst[i] = new
    
    vlst = []
    for item in eqlst:
      vlst += varlist(item)
    if len(eqlst) < len(to_compute):
        return None
    out = and0(TreeNode("f_and", eqlst), to_compute)
    if out is not None:
        return out
    if len(eqlst) == 1:
        return eqlst[0]
    return TreeNode("f_and", eqlst)
equation = None
variable = None

while cmd_mode and True:
    tmp = input(">>> ")
    try:
        orig = equation
        if tmp.split(" ")[0] == "solve":
            if out.name == "f_and":
                out = cs(equation.children, parser.take_input(tmp.split(" ")[1]))
            else:
                out = cs([equation], parser.take_input(tmp.split(" ")[1]))
            if out is not None:
                equation =  out
                print(equation)
        elif tmp.split(" ")[0] == "find":
            if out.name == "f_and":
                out = find(equation.children, [parser.take_input(item) for item in tmp.split(" ")[1:]])
            else:
                out = find([equation], [parser.take_input(item) for item in tmp.split(" ")[1:]])
            if out is not None:
                equation =  out
                print(equation)
        elif tmp.split(" ")[0] == "factor":
            if len(tmp.split(" "))==2 and tmp.split(" ")[1] == "complex":
                equation = quadratic(equation, False, None, False, True)
            elif len(tmp.split(" "))==2:
                equation = quadratic(equation, False, str_form(parser.take_input(tmp.split(" ")[1])), False)
            else:
                equation = quadratic(equation)
            equation = clearv(equation)
            equation = simplify(equation)
            equation = take_common(equation)
            print(equation)
        elif tmp == "factor1":   
            equation = quadratic(equation, True)
            equation = clearv(equation)
            equation = simplify(equation)
            print(equation)
        elif tmp == "expand":
             if equation.name == "f_eq":
                 eq2 = copy.deepcopy(equation.children[1])
                 equation = TreeNode("f_eq", [expand2(e1(equation)), eq2])
             else:
                 equation = expand2(equation)
             equation = solve(equation)
             print(equation)
        elif tmp == "simplify":
             equation = clearv(equation)
             equation = solve(equation)
             print(equation)
        elif tmp == "show":
            print(equation)
        elif tmp == "trig0":
            equation = replace_eq2(equation)
            print(equation)
        elif tmp.split(" ")[0] == "trigsolve":
            equation = trigsolve(equation, parser.take_input(tmp.split(" ")[1]))
            print(equation)
        elif tmp == "trig":
            equation = replace_eq(equation)
            equation = take_common(equation)
            equation = solve(expand2(formula_4(equation)))
            print(equation)
        elif tmp == "difsolve":
            
            tmp4 = diffsolve(equation)
            if tmp4 is None:
                tmp5 = subshomodif(equation)
                if tmp5 is None:
                    equation = lineardiff(equation)
                else:
                    equation = tmp5
            else:
                equation = tmp4
            print(equation)
        elif tmp.split(" ")[0] == "trig1":
            def command(equation):
                if equation.name == "f_mul":
                    equation = replace_eq3(equation)
                    if len(tmp.split(" "))>1 and tmp.split(" ")[1] == "allsin":
                        equation = formula_3(equation)
                    return equation
                return TreeNode(equation.name, [command(child) for child in equation.children])
            equation = solve(command(equation))
            print(equation)
        elif tmp == "apart":
            equation = decompose(equation)
            print(equation)
        elif tmp == "fraction":
            equation = fraction2(equation)
            print(equation)
        elif tmp.split(" ")[0] == "inverse":
            if equation.name == "f_eq":
                mark = parser.take_input(tmp.split(" ")[1])
                equation = TreeNode("f_eq", [inverse2(equation.children[0], mark.name), mark])
                equation = solve(expand2(equation))
            print(equation)
            equation = eqhandle(equation)
        elif tmp == "d/dx":
            equation = diffx(equation)
            print(equation)
        elif tmp[:2] == "Sd":
            plog("THOUGHT PROCESS = \n")
            var = {"x":"v_0", "y":"v_1", "z":"v_2", "c": "v_5"}[tmp[2]]
            special = TreeNode("f_int", [equation, tree_form(var)])
            su, sp, bypart = True, False, False
            if "byparts" in tmp:
                bypart = True
                su = False
            if "subs" in tmp:
                su = True
            if "graph" in tmp:
                sp = True
            n = 4
            try:
                n = int(tmp.split(" ")[-1])
            except:
                pass
            tmp2 = integratex(equation, 5, var, [], bypart, sp, su)
            tmp2 = handle_intrec(special, tmp2)
            tmp2 = solve(expand2(tmp2))
            plog(f"{'  '*tab}the solution is {str(tmp2)}\n")
            if tmp2 is None:
                print("failed to integrate")
            else:
                equation = tmp2
                print(equation)
        elif tmp.split(" ")[0] == "func":
            if variable is None:
                variable = {}
            out = normal(copy.deepcopy(tmp.split(" ")[2]), None)
            out = intreg(out)
            variable[tmp.split(" ")[1]] = out
            print(out)
        elif tmp == "mat":
            if "f_list" in str_form(equation):
                equation = mat00(equation)
            else:
                equation = and0(equation)
            print(equation)
        elif tmp == "compute":
            print(compute(equation))
        elif tmp.split(" ")[0] == "limit":
            var = parser.take_input(tmp.split(" ")[1])
            n = parser.take_input(tmp.split(" ")[2])
            eq2 = TreeNode("f_add",[tree_form("v_0"), n])
            eq = replace2(copy.deepcopy(equation), tree_form("v_0"), eq2)
            orig = equation
            tab = 0
            equation = find_limit(eq, 3, var)
            print()
            print(f"limit {var}->{n} " + str(orig) + " = " + str(equation))
        else:
            out = normal(tmp, copy.deepcopy(variable))
            equation = solve(out)
            print(equation)
    except Exception as error:
        equation = orig
        print(error)
