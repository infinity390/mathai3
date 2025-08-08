import copy
class TreeNode:
    def __init__(self, name, children=None):
        self.name = name
        self.children = children or []

    def fx(self, fxname):
        return TreeNode("f_" + fxname, [self])

    def __repr__(self):
        return string_equation(str_form(self))

    def __eq__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        elif not isinstance(other, TreeNode):
            return NotImplemented
        return str_form(self) == str_form(other)

    def __add__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return TreeNode("f_add", [self, other])

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return TreeNode("f_mul", [self, other])

    def __rmul__(self, other):
        return self.__mul__(other)

    def __sub__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return self + (tree_form("d_-1") * other)

    def __rsub__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return other + (tree_form("d_-1") * self)

    def __pow__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return TreeNode("f_pow", [self, other])

    def __rpow__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return TreeNode("f_pow", [other, self])

    def __truediv__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return self * (other ** tree_form("d_-1"))

    def __rtruediv__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return other * (self ** tree_form("d_-1"))

    def __and__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return TreeNode("f_and", [self, other])

    def __rand__(self, other):
        return self.__and__(other)

    def __or__(self, other):
        if isinstance(other, int):
            other = tree_form("d_" + str(other))
        return TreeNode("f_or", [self, other])

    def __ror__(self, other):
        return self.__or__(other)

    def __neg__(self):
        return tree_form("d_-1") * self

    def __hash__(self):
        return hash(str_form(self))

    
def str_form(node):
    def recursive_str(node, depth=0):
        result = "{}{}".format(' ' * depth, node.name)
        for child in node.children:
            result += "\n" + recursive_str(child, depth + 1)
        return result
    if not isinstance(node, TreeNode):
        return "d_"+str(node)
    return recursive_str(node)
def replace(equation, find, r):
  if str_form(equation) == str_form(find):
    return r
  col = TreeNode(equation.name, [])
  for child in equation.children:
    col.children.append(replace(child, find, r))
  return col
def tree_form(tabbed_strings):
    lines = tabbed_strings.split("\n")
    root = TreeNode("Root")
    current_level_nodes = {0: root}
    stack = [root]
    for line in lines:
        level = line.count(' ')
        node_name = line.strip()
        node = TreeNode(node_name)
        while len(stack) > level + 1:
            stack.pop()
        parent_node = stack[-1]
        parent_node.children.append(node)
        current_level_nodes[level] = node
        stack.append(node)
    return root.children[0]
def string_equation_helper(equation_tree):
    if equation_tree.children == []:
        if equation_tree.name[:2]=="g_":
            return '"'+equation_tree.name[2:]+'"'
        return equation_tree.name
    extra = ""
    
    if equation_tree.name.count("=") == 1:
        extra += "["+equation_tree.name[-1]+"]"
        
        equation_tree.name = equation_tree.name[:3]
    if equation_tree.name == "f_list":
        return "["+",".join([string_equation_helper(child) for child in equation_tree.children])+"]"
    s = "(" 
    if len(equation_tree.children) == 1 or equation_tree.name[2:] in [chr(ord("A")+i) for i in range(26)]+["int", "pdif", "dif", "A", "B", "C", "covariance"]:
        s = equation_tree.name[2:] + s
    sign = {"f_covariance": ",", "f_ge":">=", "f_le":"<=", "f_gt":">", "f_lt":"<", "f_cosec":"?" , "f_equiv": "<->", "f_sec":"?", "f_cot": "?", "f_dot": ".", "f_circumcenter":"?", "f_transpose":"?", "f_exp":"?", "f_abs":"?", "f_log":"?", "f_and":"&", "f_or":"|", "f_sub":"-", "f_neg":"?", "f_inv":"?", "f_add": "+", "f_mul": "*", "f_pow": "^", "f_poly": ",", "f_div": "/", "f_sub": "-", "f_dif": ",", "f_sin": "?", "f_cos": "?", "f_tan": "?", "f_eq": "=", "f_sqt": "?"}
    arr = []
    k = None
    if equation_tree.name not in sign.keys():
        k = ","
    else:
        k = sign[equation_tree.name]
    for child in equation_tree.children:
        arr.append(string_equation_helper(copy.deepcopy(child)))
    return s + k.join(arr) + ")"+extra
def string_equation(eq):
    
    alpha = ["x", "y", "z"]+[chr(x+ord("a")) for x in range(0,23)]
    eq = tree_form(eq)
    def rfx(eq):
        if eq.name[:2] == "f_" and eq.name[2:3] in [str(i) for i in range(26)]:
            string = ""
            if eq.name.count("=") == 1:
                string += "="+eq.name[-1]
            return TreeNode("f_"+chr(int(eq.name[2:3])+ord("A"))+string, eq.children)
        return TreeNode(eq.name, [rfx(child) for child in eq.children])
    
    eq = rfx(eq)
    for i, letter in enumerate(alpha):
        eq = replace(eq, tree_form("v_"+str(i)), tree_form(letter))
    for i in range(100, 150):
        eq = replace(eq, tree_form("v_"+str(i)), tree_form("c"+str(i-100)))
    eq = str_form(eq)
    
    eq = eq.replace("d_", "")
    eq = eq.replace("s_", "")
    eq = eq.replace("v_", "")
    eq = eq.replace("'", "")
    
    return string_equation_helper(tree_form(eq))
