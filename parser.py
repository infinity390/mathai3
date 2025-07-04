import copy
from lark import Lark, Tree
from base import *
import re
def extract_nested_lists(s, n=0):
    stack = []
    results = []
    current = ''
    indices = []
    start = None
    for i, char in enumerate(s):
        if char == '[':
            if not stack:
                start = i
            else:
                current += char
            stack.append('[')
        elif char == ']':
            stack.pop()
            if stack:
                current += char
            else:
                results.append(current)
                indices.append((start, i + 1))
                current = ''
        elif stack:
            current += char
    formatted_results = [re.sub(r'([^\[\],\s]+)', r'"\1"', lst) for lst in results]
    nested_lists = [eval(f'[{lst}]', {'__builtins__': {}}, {}) for lst in formatted_results]
    for idx, (start, end) in enumerate(reversed(indices)):
        s = s[:start] + f'{chr(len(indices)-idx-1+ord("A")+n)}' + s[end:]
    return nested_lists, s
grammar = """?start: expr
?expr: logic_equiv
?logic_equiv: logic_or
            | logic_equiv "<->" logic_or  -> equiv
?logic_or: logic_and
         | logic_or "|" logic_and  -> or
?logic_and: comparison
          | logic_and "&" comparison  -> and
?comparison: arithmetic
           | comparison "=" arithmetic  -> eq
           | comparison "<" arithmetic  -> lt
           | comparison ">" arithmetic  -> gt
           | comparison "<=" arithmetic -> le
           | comparison ">=" arithmetic -> ge
?arithmetic: term
           | arithmetic "+" term   -> add
           | arithmetic "-" term   -> sub
?term: factor
     | term "*" factor  -> mul
     | term "/" factor  -> div
     | term "." factor  -> dot
?factor: not_expr
       | factor "^" not_expr  -> pow
?not_expr: "!" not_expr  -> not
         | "-" not_expr  -> neg
         | base
?base: NUMBER              -> number
     | FUNC_NAME "(" [expr ("," expr)*] ")" index* -> func
     | VARIABLE index*     -> variable
     | "[" [expr ("," expr)*] "]" -> list
     | "(" expr ")"        -> paren
     | CNUMBER             -> cnumber
?index: "[" expr "]"     -> index
FUNC_NAME: "midpoint" | "sum" | "sum2" | "F" | "W" | "V" | "equationrhs" | "transpose" | "equationlhs" | "equation" | "error" | "covariance" | "variance" | "expect" | "A" | "B" | "C" | "mag" | "rad" | "laplace" | "diverge" | "pdif" | "gradient" | "curl" | "point1" | "point2" | "dot" | "point3" | "line1" | "line2" | "line3" | "sin" | "circumcenter" | "eqtri" | "linesegment" | "cos" | "tan" | "log" | "sqrt" | "integrate" | "dif" | "abs" | "transpose" | "cosec" | "sec" | "cot" | "arctan" | "arcsin" | "arccos" | "log10"
VARIABLE: /[a-z]/ | "nabla"
CNUMBER: /c[0-9]+/
%import common.NUMBER
%import common.WS_INLINE
%ignore WS_INLINE"""


def take_input(equation, funclist=None):
  global grammar
  equation = copy.copy(equation.replace(" ", ""))
  grammar2 = copy.deepcopy(grammar)
  if funclist is not None:
      output = grammar2.split("\n")
      for i in range(len(output)):
          if "FUNC_NAME:" in output[i]:
              output[i] = output[i].replace("FUNC_NAME: ", "FUNC_NAME: " + " | ".join(['"' + x + '"' for x in funclist]) + " | ")
      grammar2 = "\n".join(output)
  parser_main = Lark(grammar2, start='start', parser='lalr')
  parse_tree = parser_main.parse(equation)
  def convert_to_treenode(parse_tree):
      def tree_to_treenode(tree):
          if isinstance(tree, Tree):
              node = TreeNode(tree.data)
              node.children = [tree_to_treenode(child) for child in tree.children]
              return node
          else:
              return TreeNode(str(tree))
      return tree_to_treenode(parse_tree)
  def remove_past(equation):
      if equation.name in {"number", "paren", "func", "variable", "cnumber"}:
          if len(equation.children) == 1:
            for index, child in enumerate(equation.children):
              equation.children[index] = remove_past(child)
            return equation.children[0]
          else:
            for index, child in enumerate(equation.children):
              equation.children[index] = remove_past(child)
            return TreeNode(equation.children[0].name, equation.children[1:])
      coll = TreeNode(equation.name, [])
      for child in equation.children:
          coll.children.append(remove_past(child))
      return coll
  def prefixindex(equation):
      out = None
      for index, child in enumerate(equation.children):
          if child.name == "index":
              out = str(child.children[0])
              equation.children.pop(index)
              break
      if out is not None:
          return TreeNode(equation.name + "="+out, equation.children)
      return TreeNode(equation.name, [prefixindex(child) for child in equation.children])
  tree_node = convert_to_treenode(parse_tree)
  tree_node = remove_past(tree_node)
  tree_node = prefixindex(tree_node)
  
  def fxchange(tree_node):
    nonlocal funclist
    tmp3 = []
    if funclist is not None:
        tmp3 = funclist
    return TreeNode("f_"+tree_node.name if tree_node.name in tmp3+["sum", "sum2", "V", "W", "transpose", "equationrhs", "equationlhs", "equation", "covariance", "variance", "expect", "error", "laplace", "dot", "curl", "pdif", "diverge", "gradient", "rad", "ge", "le", "gt", "lt", "F", "rad", "eqtri", "linesegment", "midpoint", "mag", "point1", "point2", "point3", "line1", "line2", "line3", "log10", "arcsin", "arccos", "arctan", "list", "cosec", "sec", "cot", "equiv", "or", "not", "and", "circumcenter", "transpose", "eq", "sub", "neg", "inv", "add", "sin", "cos", "tan", "mul", "integrate", "dif", "pow", "div", "log", "abs"]\
                    else "d_"+tree_node.name, [fxchange(child) for child in tree_node.children])
  tree_node = fxchange(tree_node)
  tree_node = replace(tree_node, tree_form("d_e"), tree_form("s_e"))
  tree_node = replace(tree_node, tree_form("d_pi"), tree_form("s_pi"))
  tree_node = replace(tree_node, tree_form("d_i"), tree_form("s_i"))
  tree_node = replace(tree_node, tree_form("d_nabla"), tree_form("s_nabla"))

  def rfx(tree_node):
      if tree_node.name[:2] == "d_" and tree_node.name[2:3] in ["A", "B", "C"]:
          out = "f_"+str(ord(tree_node.name[2:3])-ord("A"))
          if tree_node.name.count("=") == 1:
              out += "=" + tree_node.name[-1]
          return TreeNode(out, tree_node.children)
      if tree_node.name[:2] == "d_" and tree_node.name[2:3] in ["a", "b", "c"]:
          out = "v_"+str(ord(tree_node.name[2:3])-ord("a"))
          if tree_node.name.count("=") == 1:
              out += "=" + tree_node.name[-1]
          return TreeNode(out, tree_node.children)
      if tree_node.name[:3] == "d_c":
          return tree_form("v_" + str(int(tree_node.name[3:])+100))
      return TreeNode(tree_node.name, [rfx(child) for child in tree_node.children])

  for i in range(26):
    alpha = ["x", "y", "z"]+[chr(x+ord("a")) for x in range(0,23)]
    tree_node = replace(tree_node, tree_form("d_"+alpha[i]), tree_form("v_"+str(i)))
  tmp5 = rfx(tree_node)
  return tmp5
