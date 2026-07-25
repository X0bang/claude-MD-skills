#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MD 建模前结构体检 / Pre-MD structure sanity check

零依赖(纯标准库),任何有 python3 的机器都能跑。
Zero dependencies (stdlib only) — runs anywhere python3 exists.

分级原则 / Tier principle:
  A 级 = 跑 MD 永远修不好的错误(手性、顺式肽键、缺原子、穿环)
         必须在建模前修掉,否则错误会一路带进所有下游结果。
  B 级 = 能量极小化能修好的问题(键长、重叠、断链)
         会让建库工具报出天文数字的单点能量,但不影响最终结果。
  C 级 = 文件格式问题,会让上传/解析直接失败。
  D 级 = 膜蛋白专属(拓扑朝向、疏水失配、positive-inside)。
  E 级 = 多体系横向一致性 —— 单看每个都正常、放一起才发现不对的那类问题。

用法 / Usage:
  check_structure.py <file.pdb> [more.pdb ...]
  check_structure.py *.pdb --membrane            # 打开 D 级膜蛋白检查
  check_structure.py *.pdb --membrane --normal z # 指定膜法向(默认 z)
  check_structure.py *.pdb --quiet               # 只报问题,不报通过项

退出码 / Exit codes: 0=无 A/C 级问题  1=有 A 级  2=有 C 级  3=两者都有
"""

import sys, os, math, argparse
from collections import OrderedDict, defaultdict

# ─────────────────────────────────────────────────────────────────────────────
# 参考数据 / Reference data
# ─────────────────────────────────────────────────────────────────────────────

# 标准氨基酸重原子 + 侧链连接性(主链 N-CA / CA-C / C-O 所有残基共有)
TEMPLATES = {
    "ALA": (["CB"], [("CA", "CB")]),
    "ARG": (["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"),
             ("CZ", "NH1"), ("CZ", "NH2")]),
    "ASN": (["CB", "CG", "OD1", "ND2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")]),
    "ASP": (["CB", "CG", "OD1", "OD2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")]),
    "CYS": (["CB", "SG"], [("CA", "CB"), ("CB", "SG")]),
    "GLN": (["CB", "CG", "CD", "OE1", "NE2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")]),
    "GLU": (["CB", "CG", "CD", "OE1", "OE2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")]),
    "GLY": ([], []),
    "HIS": (["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"),
             ("ND1", "CE1"), ("CD2", "NE2"), ("CE1", "NE2")]),
    "ILE": (["CB", "CG1", "CG2", "CD1"],
            [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")]),
    "LEU": (["CB", "CG", "CD1", "CD2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")]),
    "LYS": (["CB", "CG", "CD", "CE", "NZ"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")]),
    "MET": (["CB", "CG", "SD", "CE"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")]),
    "PHE": (["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
             ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")]),
    "PRO": (["CB", "CG", "CD"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N")]),
    "SER": (["CB", "OG"], [("CA", "CB"), ("CB", "OG")]),
    "THR": (["CB", "OG1", "CG2"], [("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")]),
    "TRP": (["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "NE1"),
             ("NE1", "CE2"), ("CD2", "CE2"), ("CD2", "CE3"), ("CE2", "CZ2"),
             ("CE3", "CZ3"), ("CZ2", "CH2"), ("CZ3", "CH2")]),
    "TYR": (["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
            [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"),
             ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CZ", "OH")]),
    "VAL": (["CB", "CG1", "CG2"], [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")]),
}
BACKBONE = ["N", "CA", "C", "O"]
BB_BONDS = [("N", "CA"), ("CA", "C"), ("C", "O")]

# 端基封端(patch)的成键。MD 体系几乎都封端,不认这些会把封端原子误报成"原子重叠"。
#   ACE 乙酰封端 (CHARMM): CAY-CY-OY,CY 连主链 N
#   CT2/CT3 酰胺封端 (CHARMM): C-NT,CT3 另有 NT-CAT
#   标准羧基端: C-OXT
PATCH_BONDS = [("CAY", "CY"), ("CY", "OY"), ("CY", "N"),
               ("C", "NT"), ("NT", "CAT"), ("C", "OXT"),
               ("CH3", "C"), ("N", "CH3")]      # AMBER: ACE 的 CH3-C,NME 的 N-CH3
PATCH_IDEAL = {frozenset(("CAY", "CY")): 1.52, frozenset(("CY", "OY")): 1.23,
               frozenset(("CY", "N")): 1.35, frozenset(("C", "NT")): 1.33,
               frozenset(("NT", "CAT")): 1.45, frozenset(("C", "OXT")): 1.25,
               frozenset(("CH3", "C")): 1.52}

# CHARMM / AMBER 常见残基名别名 -> 标准名
ALIASES = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS",
           "HIP": "HIS", "CYX": "CYS", "CYM": "CYS", "ASH": "ASP", "GLH": "GLU",
           "LYN": "LYS", "MSE": "MET"}

# 原子名别名:CHARMM 与 PDB 标准的差异。不归一化会把"命名不同"误报成"缺原子"。
#   ILE 的 δ 碳:CHARMM 叫 CD,PDB 标准叫 CD1
#   C 端羧基氧:CHARMM 叫 OT1/OT2,PDB 标准叫 O/OXT
#   硒代甲硫氨酸 MSE:硫位上是硒,原子名 SE(归一化成 MET 后模板要的是 SD)
ATOM_ALIASES = {
    "ILE": {"CD": "CD1"},
    "MET": {"SE": "SD"},          # MSE -> MET 之后 SE 要映射成 SD
    "*":   {"OT1": "O", "OT2": "OXT", "O1": "O", "O2": "OXT"},
}


def norm_atom(resname, atom):
    rn = ALIASES.get(resname, resname)
    return ATOM_ALIASES.get(rn, {}).get(atom, ATOM_ALIASES["*"].get(atom, atom))

# 侧链净电荷。必须用**原始残基名**查表 —— 先经 std() 归一化会把 HSP/ASH/LYN 等
# 质子化态信息抹掉,导致电荷算错(HSP 会变成 HIS 而得到 0 而非 +1)。
CHARGE = {"ARG": +1, "LYS": +1, "ASP": -1, "GLU": -1,
          "HSP": +1, "HIP": +1,                      # 双质子化 His
          "ASH": 0, "GLH": 0, "LYN": 0, "CYM": -1}   # 非标准质子化态


def res_charge(raw_name):
    if raw_name in CHARGE:
        return CHARGE[raw_name]
    return CHARGE.get(std(raw_name), 0)

# 理想键长 (Å) 与容差
IDEAL_BONDS = {("N", "CA"): 1.458, ("CA", "C"): 1.525, ("C", "O"): 1.231}
IDEAL_SIDECHAIN = 1.52     # 通用 sp3 C-C
BOND_TOL = 0.15

# 侧链具体键长 —— 用通用 1.52 会把芳香环(1.39)和胍基/酰胺 C-N(1.33)全部误报。
IDEAL_SC = {}
def _sc(res, pairs, val):
    for a, b in pairs:
        IDEAL_SC[(res, frozenset((a, b)))] = val

for _r in ("PHE", "TYR"):
    _sc(_r, [("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"),
             ("CE1", "CZ"), ("CE2", "CZ")], 1.39)
_sc("TYR", [("CZ", "OH")], 1.38)
_sc("TRP", [("CG", "CD1")], 1.37); _sc("TRP", [("CG", "CD2")], 1.43)
_sc("TRP", [("CD1", "NE1"), ("NE1", "CE2")], 1.37)
_sc("TRP", [("CD2", "CE2"), ("CD2", "CE3"), ("CE2", "CZ2"),
            ("CE3", "CZ3"), ("CZ2", "CH2"), ("CZ3", "CH2")], 1.40)
_sc("HIS", [("CG", "ND1")], 1.38); _sc("HIS", [("CG", "CD2")], 1.35)
_sc("HIS", [("ND1", "CE1"), ("CE1", "NE2")], 1.32); _sc("HIS", [("CD2", "NE2")], 1.37)
_sc("ARG", [("CD", "NE")], 1.46); _sc("ARG", [("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")], 1.33)
_sc("LYS", [("CE", "NZ")], 1.49)
_sc("ASN", [("CG", "OD1")], 1.23); _sc("ASN", [("CG", "ND2")], 1.33)
_sc("GLN", [("CD", "OE1")], 1.23); _sc("GLN", [("CD", "NE2")], 1.33)
_sc("ASP", [("CG", "OD1"), ("CG", "OD2")], 1.25)
_sc("GLU", [("CD", "OE1"), ("CD", "OE2")], 1.25)
_sc("SER", [("CB", "OG")], 1.42); _sc("THR", [("CB", "OG1")], 1.43)
_sc("CYS", [("CB", "SG")], 1.81)
_sc("MET", [("CG", "SD"), ("SD", "CE")], 1.81)
PEPTIDE_CN = 1.329
CA_CA_TRANS = 3.80

# 芳香环(用于平面性与穿环检测)
RINGS = {
    "PHE": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
    "TYR": [["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]],
    "HIS": [["CG", "ND1", "CE1", "NE2", "CD2"]],
    "TRP": [["CG", "CD1", "NE1", "CE2", "CD2"],
            ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]],
    "PRO": [["N", "CA", "CB", "CG", "CD"]],
}

# Kyte-Doolittle 疏水性(用于识别跨膜段)
KD = {"ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
      "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3,
      "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
      "LYS": -3.9, "ARG": -4.5}

# ─────────────────────────────────────────────────────────────────────────────
# 纯 Python 向量运算 / Pure-python vector math
# ─────────────────────────────────────────────────────────────────────────────

def sub(a, b):   return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def add(a, b):   return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def scale(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def dot(a, b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a):     return math.sqrt(dot(a, a))
def dist(a, b):  return norm(sub(a, b))

def unit(a):
    n = norm(a)
    return (0.0, 0.0, 0.0) if n < 1e-12 else scale(a, 1.0/n)

def det3(a, b, c):
    """行列式 = a · (b × c),用于手性判定"""
    return dot(a, cross(b, c))

def dihedral(p0, p1, p2, p3):
    """二面角(度),范围 -180..180"""
    b0, b1, b2 = sub(p0, p1), sub(p2, p1), sub(p3, p2)
    b1u = unit(b1)
    v = sub(b0, scale(b1u, dot(b0, b1u)))
    w = sub(b2, scale(b1u, dot(b2, b1u)))
    x, y = dot(v, w), dot(cross(b1u, v), w)
    return math.degrees(math.atan2(y, x))

def _best_normal(pts, c=None):
    """最佳拟合平面的法向 = 协方差矩阵最小特征值对应的特征向量。

    直接对 (trace(M)·I − M) 做幂迭代:它的最大特征向量正是 M 的最小特征向量。
    这样避开了"先求最大、再正交去除两次"的做法 —— 后者在**完全对称的共面点集**
    (理想正 n 边形)上会遇到 λ1 = λ2 的简并,去除后残量精确为零,只能退回一个
    硬编码方向,算出来的"离面量"其实是环的面内展宽,导致平面性误报。
    """
    n = len(pts)
    if c is None:
        c = scale((sum(p[0] for p in pts), sum(p[1] for p in pts),
                   sum(p[2] for p in pts)), 1.0/n)
    q = [sub(p, c) for p in pts]
    M = [[sum(v[i]*v[j] for v in q) for j in range(3)] for i in range(3)]
    tr = M[0][0] + M[1][1] + M[2][2]
    N = [[(tr if i == j else 0.0) - M[i][j] for j in range(3)] for i in range(3)]
    v = (0.5773, 0.5774, 0.5775)
    for _ in range(200):
        w = tuple(sum(N[i][j]*v[j] for j in range(3)) for i in range(3))
        if norm(w) < 1e-15:
            break
        v = unit(w)
    return v


def plane_rmsd(pts):
    """点集到最佳拟合平面的 RMSD"""
    n = len(pts)
    c = scale((sum(p[0] for p in pts), sum(p[1] for p in pts), sum(p[2] for p in pts)), 1.0/n)
    nrm = _best_normal(pts, c)
    return math.sqrt(sum(dot(sub(p, c), nrm)**2 for p in pts) / n)

def principal_axis(pts):
    """最大惯量主轴(幂迭代)"""
    n = len(pts)
    c = scale((sum(p[0] for p in pts), sum(p[1] for p in pts), sum(p[2] for p in pts)), 1.0/n)
    q = [sub(p, c) for p in pts]
    M = [[sum(v[i]*v[j] for v in q) for j in range(3)] for i in range(3)]
    v = (0.5773, 0.5774, 0.5775)
    for _ in range(200):
        v = unit(tuple(sum(M[i][j]*v[j] for j in range(3)) for i in range(3)))
    return v, c

# ─────────────────────────────────────────────────────────────────────────────
# PDB / GRO 解析 —— 保留原始行,以便做格式检查
# ─────────────────────────────────────────────────────────────────────────────

class Atom(object):
    __slots__ = ("serial", "name", "altloc", "resname", "chain", "resseq", "icode",
                 "xyz", "occ", "bfac", "element", "segid", "lineno", "raw")

class Structure(object):
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        # 多个不同目录下同名文件(如四个 step5_input.pdb)在 E 级表里必须能区分
        par = os.path.basename(os.path.dirname(os.path.abspath(path)))
        self.label = "%s/%s" % (par, self.name) if par else self.name
        self.atoms = []
        self.raw_lines = []
        self.fmt = "gro" if path.endswith(".gro") else "pdb"
        self.has_hydrogen = False
        self.extra_models = False        # 是否丢弃了额外的 NMR 模型
        self.bad_coord_lines = []        # 坐标列无法解析的行号
        self.dropped_altloc = 0          # 被过滤掉的次要构象原子数
        self.lost_atoms = 0              # 因链标识冲突被静默覆盖掉的原子数
        self.gro_split = 0               # .gro 按几何自动切出的链数
        (self._read_gro if self.fmt == "gro" else self._read_pdb)()
        self.residues = self._group()

    def _read_pdb(self):
        model = 0
        with open(self.path) as fh:
            for i, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                self.raw_lines.append((i, line))
                # 多模型(NMR ensemble):只取第一个模型。不处理的话后面的模型会
                # 覆盖前面的坐标,并制造满屏假"重复原子"。
                if line.startswith("MODEL"):
                    model += 1
                    if model > 1:
                        self.extra_models = True
                    continue
                if model > 1:
                    continue
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                a = Atom()
                a.lineno, a.raw = i, line
                pad = line.ljust(80)
                try:
                    a.serial = int(pad[6:11])
                except ValueError:
                    a.serial = len(self.atoms) + 1
                a.name    = pad[12:16].strip()
                a.altloc  = pad[16]
                a.resname = pad[17:20].strip()
                a.chain   = pad[21]
                try:
                    a.resseq = int(pad[22:26])
                except ValueError:
                    a.resseq = -9999
                a.icode = pad[26]
                try:
                    a.xyz = (float(pad[30:38]), float(pad[38:46]), float(pad[46:54]))
                except ValueError:
                    # 坐标列被截断或溢出(如 -1234.567 占 9 列)。记下行号继续,
                    # 而不是让整个文件解析失败。
                    self.bad_coord_lines.append(i)
                    continue
                a.occ  = _f(pad[54:60], 1.0)
                a.bfac = _f(pad[60:66], 0.0)
                a.segid   = pad[72:76].strip()
                a.element = pad[76:78].strip()
                if a.element == "H" or (not a.element and a.name[:1] == "H"):
                    self.has_hydrogen = True
                self.atoms.append(a)

    def _read_gro(self):
        with open(self.path) as fh:
            lines = fh.read().splitlines()
        n = int(lines[1])
        for i, line in enumerate(lines[2:2+n], 3):
            self.raw_lines.append((i, line))
            a = Atom()
            a.lineno, a.raw = i, line
            a.resseq  = int(line[0:5])
            a.resname = line[5:10].strip()
            a.name    = line[10:15].strip()
            a.serial  = int(line[15:20])
            # gro 单位是 nm -> Å
            a.xyz = tuple(float(line[20+8*k:28+8*k]) * 10.0 for k in range(3))
            a.altloc, a.chain, a.icode = " ", " ", " "
            a.occ, a.bfac, a.element, a.segid = 1.0, 0.0, "", ""
            if a.name[:1] == "H":
                self.has_hydrogen = True
            self.atoms.append(a)

    def _group(self):
        """按 (链标识, resseq, icode) 分组;跳过氢和水/离子/脂质

        ⚠️ 链标识必须同时看 chain ID 和 segid。CHARMM / CHARMM-GUI / psfgen 风格的
           文件里,多条链的 chain ID 可能全都一样(甚至为空),真正区分它们的是
           73-76 列的 segid。只按 chain ID 分组会把多个原体**静默合并**成一条链 ——
           后写入的覆盖先写入的,大部分结构凭空消失,而且会伪造出满屏"重复原子"。
        """
        SKIP = {"HOH", "WAT", "TIP3", "SOL", "NA", "CL", "K", "POT", "CLA", "SOD",
                "POPC", "POPE", "POPS", "POPG", "CHL1", "DPPC"}
        def is_h(a):
            if a.element in ("H", "D"):
                return True
            if not a.element and a.name[:1] in ("H", "D"):
                return True
            return a.name[:1] in ("H", "D")

        res = OrderedDict()
        heavy_in = 0             # 应当被分组的重原子数,用于事后校验有没有静默丢失
        prev = None
        gro_chain = 0
        # 先扫一遍决定每个残基采用哪个 altLoc:优先主构象,若某残基**整个**只以
        # 次要构象存在(B/C…),就回退到字母最小的那个 —— 否则该残基会静默消失,
        # 表面上只看到一处莫名其妙的"主链断裂"。
        alt_of = {}
        for a in self.atoms:
            if a.resname in SKIP or is_h(a):
                continue
            k = (chain_key(a), a.resseq, a.icode)
            cur = alt_of.get(k)
            al = a.altloc if a.altloc.strip() else " "
            if cur is None or (cur != " " and (al == " " or al < cur)):
                alt_of[k] = al
        for a in self.atoms:
            # .gro 没有 chain ID 也没有 segid。多链体系里残基号会重新从 1 开始,
            # 用这个跳变切分链 —— 否则多条链会被静默并成一条(残基号重复时甚至
            # 会互相覆盖,大部分原子凭空消失)。
            if self.fmt == "gro":
                if prev is not None and a.resseq < prev:
                    gro_chain += 1
                prev = a.resseq
                a.segid = "seg%d" % gro_chain if gro_chain else ""
            if a.resname in SKIP:
                continue
            # 氘(D)也要滤掉:中子衍射结构用元素符号 D,只按 H 过滤会把它们当成
            # 重原子,制造成百上千的假"原子重叠"。
            if is_h(a):
                continue
            key = (chain_key(a), a.resseq, a.icode)
            # 交替构象:每个残基只取选定的那一套,否则两套构象被并成一个残基,
            # 制造大量假重叠(甚至残基种类都可能不同,如 LEU/ILE 微观异质性)。
            al = a.altloc if a.altloc.strip() else " "
            if al != alt_of.get(key, " "):
                self.dropped_altloc += 1
                continue
            heavy_in += 1
            res.setdefault(key, {"name": a.resname, "atoms": OrderedDict()})
            res[key]["atoms"][norm_atom(a.resname, a.name)] = a.xyz
        # 不变量:进入分组的重原子数必须等于分组后的原子总数。不相等说明有原子
        # 被同键覆盖 —— 典型场景是多条链共用链标识且残基号重复,后写的盖掉先写的。
        grouped = sum(len(r["atoms"]) for r in res.values())
        if grouped < heavy_in:
            self.lost_atoms = heavy_in - grouped

        # .gro 完全没有链信息,而且 CHARMM-GUI 的 gro 是**跨链连续编号**的,
        # 所以残基号跳变这一招也失效。只能用几何:CA-CA 超过 4.2 Å 就是新链起点。
        # 不切的话多个原体会被当成一条链 —— 表现为假的"主链断裂",而且 D 级会
        # 对着一条拼接出来的假链算主轴和拓扑,结论完全没有意义。
        if self.fmt == "gro" and res:
            out = OrderedDict()
            ci, prev_ca = 0, None
            for (ck, rs, ic), r in res.items():
                ca = r["atoms"].get("CA")
                if prev_ca is not None and ca is not None and dist(prev_ca, ca) > 4.2:
                    ci += 1
                    self.gro_split = ci + 1
                if ca is not None:
                    prev_ca = ca
                out[("seg%d" % ci, rs, ic)] = r
            res = out
        return res

    def chains(self):
        out = OrderedDict()
        for (ch, rs, ic), r in self.residues.items():
            out.setdefault(ch, []).append(((ch, rs, ic), r))
        return out


def _f(s, d):
    try:
        return float(s)
    except ValueError:
        return d


def std(resname):
    return ALIASES.get(resname, resname)


def chain_key(a):
    """链标识 = chain ID + segid。两处(分组、重复原子检测)必须用同一套判据,
    否则 CHARMM 风格文件里靠 segid 区分的多条链会被当成同一条。"""
    ck = a.chain.strip()
    if a.segid:
        return "%s:%s" % (ck, a.segid) if ck else a.segid
    return ck


# ─────────────────────────────────────────────────────────────────────────────
# 检查项 / Checks
# ─────────────────────────────────────────────────────────────────────────────

class Report(object):
    def __init__(self, name):
        self.name = name
        self.items = defaultdict(list)      # tier -> [(status, title, detail)]

    def add(self, tier, status, title, detail=""):
        self.items[tier].append((status, title, detail))

    def counts(self, tier):
        bad = sum(1 for s, _, _ in self.items[tier] if s == "FAIL")
        warn = sum(1 for s, _, _ in self.items[tier] if s == "WARN")
        return bad, warn


# ── A 级 ─────────────────────────────────────────────────────────────────────

def check_chirality(st, rep):
    """A1/A2 手性 —— MD 永远修不好

    判据:
      CA:      det(N-CA,  C-CA,   CB-CA)  > 0  => L 型
      ILE-CB:  det(CA-CB, CG1-CB, CG2-CB) > 0  => 正确构型 (2S,3S)
      THR-CB:  det(CA-CB, OG1-CB, CG2-CB) > 0  => 正确构型 (2S,3R)

    符号依据(三条独立证据,已交叉验证):
      1. CIP 推导:对优先级 p1>p2>p3>p4 的手性中心,S 型 <=> det[p1,p2,p3] > 0。
         CA 处优先级 N > C' > CB > HA,L-氨基酸为 (S),故 det > 0 = L。
         ⚠️ Cys 是个坑:S(16) > O(8) 使 CB > C',CIP 字母翻成 (R) —— 但那只是
         命名翻转,把行列式的后两列换回同样的原子顺序后符号不变,20 种残基统一适用。
      2. 实验结构:1UBQ(1.8 Å)/1CRN(1.5 Å)/4LZT(0.95 Å) 共 229 个带 CB 的残基,
         行列式 100% 为正(均值 +2.5,最小 +1.93),其中 14 个 CYS 同样全正。
      3. CHARMM36 拓扑 top_all36_prot.rtf 的 IC 表:L-Ala 的 improper
         N-C-CA-CB = +123.23°,据此重建坐标得 det = +2.467,镜像得 -2.467。

    ⚠️ 不要用 "CHARMM improper 能量接近 0" 来论证手性正确 —— CHARMM36 蛋白力场
       根本没有 CA 手性 improper(只有 IMPR N -C CA HN 和 IMPR C CA +N O,约束的是
       酰胺 H 与羰基 O 的平面性)。D 型残基给出的 improper 能量完全相同。

    判别余量很大:真实 L 残基的归一化行列式在 +0.58~+0.83,离阈值 0 极远,无边界误判风险。
    """
    d_ca, d_cb = [], []
    for (ch, rs, ic), r in st.residues.items():
        rn, at = std(r["name"]), r["atoms"]
        if all(k in at for k in ("N", "CA", "C", "CB")):
            v = det3(sub(at["N"], at["CA"]), sub(at["C"], at["CA"]), sub(at["CB"], at["CA"]))
            if v < 0:
                d_ca.append("%s%d%s" % (ch.strip() or "-", rs, rn))
        if rn == "ILE" and all(k in at for k in ("CA", "CB", "CG1", "CG2")):
            if det3(sub(at["CA"], at["CB"]), sub(at["CG1"], at["CB"]), sub(at["CG2"], at["CB"])) < 0:
                d_cb.append("%s%d ILE-CB" % (ch.strip() or "-", rs))
        if rn == "THR" and all(k in at for k in ("CA", "CB", "OG1", "CG2")):
            if det3(sub(at["CA"], at["CB"]), sub(at["OG1"], at["CB"]), sub(at["CG2"], at["CB"])) < 0:
                d_cb.append("%s%d THR-CB" % (ch.strip() or "-", rs))

    if d_ca:
        rep.add("A", "FAIL", "手性错误:%d 个 D 型氨基酸 (CA 中心)" % len(d_ca),
                "位置: %s\n    自然界蛋白质只有 L 型。补氢时会把氢放进已被占据的位置,\n"
                "    导致单点能量爆到千万级。跑 MD 改不了手性(要断键重接),必须建模前修。"
                % ", ".join(d_ca))
    else:
        rep.add("A", "PASS", "手性 (CA):全部 L 型")

    if d_cb:
        rep.add("A", "FAIL", "手性错误:%d 个 Ile/Thr 第二手性中心 (CB)" % len(d_cb),
                "位置: %s\n    Ile 和 Thr 的 CB 也是手性中心,同样 MD 修不好。" % ", ".join(d_cb))
    elif any(std(r["name"]) in ("ILE", "THR") for r in st.residues.values()):
        rep.add("A", "PASS", "手性 (Ile/Thr 的 CB):全部正确")


def check_cis_peptide(st, rep):
    """A3 顺式肽键 —— 非脯氨酸的顺式肽键几乎一定是建模错误,MD 翻不过来"""
    cis, cis_pro = [], []
    for ch, lst in st.chains().items():
        for i in range(len(lst) - 1):
            (_, rs0, _), r0 = lst[i]
            (_, rs1, _), r1 = lst[i + 1]
            a0, a1 = r0["atoms"], r1["atoms"]
            if not all(k in a0 for k in ("CA", "C")) or not all(k in a1 for k in ("N", "CA")):
                continue
            if dist(a0["C"], a1["N"]) > 2.0:       # 断链,不是肽键
                continue
            omega = dihedral(a0["CA"], a0["C"], a1["N"], a1["CA"])
            if abs(omega) < 90:
                tag = "%s%d-%d (ω=%.0f°)" % (ch.strip() or "-", rs0, rs1, omega)
                (cis_pro if std(r1["name"]) == "PRO" else cis).append(tag)
    if cis:
        rep.add("A", "FAIL", "顺式肽键(非脯氨酸)%d 处" % len(cis),
                "位置: %s\n    正常肽键是反式(ω≈180°)。非 Pro 的顺式极罕见,基本是建模错误。\n"
                "    MD 翻转肽键的能垒极高,常规模拟时间内翻不过来。" % ", ".join(cis))
    if cis_pro:
        rep.add("A", "WARN", "顺式脯氨酸 %d 处" % len(cis_pro),
                "位置: %s\n    Pro 的顺式肽键天然存在(约 5%%),如果是实验结构则正常;\n"
                "    如果是建模产物,建议确认是否有依据。" % ", ".join(cis_pro))
    if not cis and not cis_pro:
        rep.add("A", "PASS", "肽键构型:全部反式")


def check_missing_atoms(st, rep):
    """A4 缺失重原子"""
    missing, unknown = [], []
    for (ch, rs, ic), r in st.residues.items():
        rn = std(r["name"])
        if rn not in TEMPLATES:
            unknown.append("%s%d %s" % (ch.strip() or "-", rs, r["name"]))
            continue
        want = set(BACKBONE) | set(TEMPLATES[rn][0])
        have = set(r["atoms"])
        lack = want - have
        # 注意:不要无条件丢弃 "O"。C 端的 OT1/O1 已由 ATOM_ALIASES 归一化成 O,
        # 所以这里缺 O 就是真的缺主链羰基氧 —— 那是 A 级问题(补出来的 O 方向由
        # φ/ψ 决定,可能完全错),必须报出来。
        if lack:
            missing.append("%s%d%s 缺 %s" % (ch.strip() or "-", rs, rn, ",".join(sorted(lack))))
    if missing:
        rep.add("A", "FAIL", "缺失重原子:%d 个残基" % len(missing),
                "\n    ".join(missing[:15]) + ("\n    ..." if len(missing) > 15 else "") +
                "\n    建库工具通常会自动补,但补出来的位置未必合理;若缺主链原子必须先修。")
    else:
        rep.add("A", "PASS", "重原子完整")
    if unknown:
        rep.add("A", "WARN", "非标准残基 %d 个" % len(unknown),
                ", ".join(unknown[:15]) + "\n    需要确认力场是否有对应参数。")


def check_planarity(st, rep):
    """A5 芳香环平面性 —— 环翘了是建模错误,力场也拉不回来"""
    bad = []
    for (ch, rs, ic), r in st.residues.items():
        rn = std(r["name"])
        if rn not in RINGS or rn == "PRO":
            continue
        for ring in RINGS[rn]:
            pts = [r["atoms"][n] for n in ring if n in r["atoms"]]
            if len(pts) < len(ring):
                continue
            d = plane_rmsd(pts)
            if d > 0.10:
                bad.append("%s%d%s (偏离 %.2f Å)" % (ch.strip() or "-", rs, rn, d))
    if bad:
        rep.add("A", "FAIL", "芳香环不平面:%d 处" % len(bad),
                ", ".join(bad) + "\n    芳香环必须共平面(RMSD < 0.1 Å)。翘曲说明建模时几何没约束好。")
    else:
        rep.add("A", "PASS", "芳香环平面性正常")


def check_ring_piercing(st, rep):
    """A6 侧链穿环 —— 键穿过环内部,拓扑上打了结,极小化只会把它锁死"""
    rings = []
    for (ch, rs, ic), r in st.residues.items():
        rn = std(r["name"])
        for ring in RINGS.get(rn, []):
            pts = [r["atoms"][n] for n in ring if n in r["atoms"]]
            if len(pts) == len(ring):
                c = scale((sum(p[0] for p in pts), sum(p[1] for p in pts),
                           sum(p[2] for p in pts)), 1.0/len(pts))
                # 法向用最佳拟合平面,不用三点叉积 —— 三点法会给环上原子随机的
                # 微小离面量,让"两端异侧"频繁假成立。
                nv = _best_normal(pts, c)
                # "在环内"的判据要用内切半径(apothem),不是外接半径。
                # 正 n 边形 apothem/外接 = cos(π/n):六元环 0.866、五元环 0.809。
                # 原先硬编码的 0.85 大于五元环的 0.809,任何五元环的边都会被误判。
                mean_rad = sum(dist(p, c) for p in pts) / len(pts)
                rad = mean_rad * math.cos(math.pi / len(ring))
                rings.append((ch, rs, rn, c, nv, rad, set(id(p) for p in pts), pts))

    bonds = _bond_list(st)
    hits = []
    for (ch, rs, rn, c, nv, rad, ids, pts) in rings:
        for (p, q, tag, _bt) in bonds:
            if dist(c, p) > rad + 3.0 and dist(c, q) > rad + 3.0:
                continue
            # 环自身的边不算穿环。TRP 是并环(吲哚),五元环与六元环共用 CD2-CE2 —— 这条边
            # 落在环平面内,而五元环边中点到环心约 0.81*rad,正好在下面 0.85*rad 的阈值内,
            # 浮点噪声就会把它误判成穿环。
            if id(p) in ids or id(q) in ids:
                continue
            dp, dq = dot(sub(p, c), nv), dot(sub(q, c), nv)
            if dp * dq >= 0:                    # 两端同侧,没穿过平面
                continue
            # 注意:不要再加"近共面就跳过"的过滤。排除环上原子(上面那步)已经足够消除
            # 全部误报(在 1UBQ/1CRN/3NIR/1A2P/4HHB/5PTI/2AXT 共 36 条链 5396 个残基上
            # 零误报);而 |dp| 阈值会让**倾角小于约 30° 的浅穿环**静默漏检 —— 而环肽穿环、
            # 脂质尾巴穿 Trp 恰恰多是浅角度的。
            t = dp / (dp - dq)
            x = add(p, scale(sub(q, p), t))
            if dist(x, c) < rad:                # rad 已是内切半径
                hits.append("%s 环 %s%d%s 被 %s 穿过" % ("", ch.strip() or "-", rs, rn, tag))
    if hits:
        rep.add("A", "FAIL", "穿环 %d 处" % len(hits),
                "\n    ".join(sorted(set(hits))[:10]) +
                "\n    化学键穿过了环的内部,拓扑上无法解开。极小化只会把这个结锁得更死。")
    else:
        rep.add("A", "PASS", "无穿环")


def check_disulfide(st, rep):
    """A7 二硫键几何"""
    sgs = [((ch, rs), r["atoms"]["SG"]) for (ch, rs, ic), r in st.residues.items()
           if std(r["name"]) == "CYS" and "SG" in r["atoms"]]
    if not sgs:
        return
    close, odd = [], []
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            d = dist(sgs[i][1], sgs[j][1])
            if d < 2.5:
                close.append("%s%d-%s%d (%.2f Å)" % (sgs[i][0][0].strip() or "-", sgs[i][0][1],
                                                     sgs[j][0][0].strip() or "-", sgs[j][0][1], d))
            elif d < 3.2:
                odd.append("%s%d-%s%d (%.2f Å)" % (sgs[i][0][0].strip() or "-", sgs[i][0][1],
                                                   sgs[j][0][0].strip() or "-", sgs[j][0][1], d))
    if close:
        rep.add("A", "PASS", "检出二硫键 %d 对" % len(close), ", ".join(close) +
                "\n    建库时必须显式声明,否则会被当成两个自由巯基。")
    if odd:
        rep.add("A", "WARN", "Cys 距离暧昧 %d 对" % len(odd),
                ", ".join(odd) + "\n    在 2.5–3.2 Å 之间,既不像成键也不像分开,需人工判断。")


# ── B 级 ─────────────────────────────────────────────────────────────────────

def _bond_list(st):
    """按残基模板生成键表(不依赖距离,所以坏结构也能正确判定)"""
    out = []
    for ch, lst in st.chains().items():
        for i, ((_, rs, ic), r) in enumerate(lst):
            rn, at = std(r["name"]), r["atoms"]
            if rn not in TEMPLATES:
                continue          # 配体/非标准残基没有模板,套主链理想值只会误报
            pairs = list(BB_BONDS) + list(TEMPLATES.get(rn, ([], []))[1]) + PATCH_BONDS
            for a, b in pairs:
                if a in at and b in at:
                    out.append((at[a], at[b],
                                "%s%d%s %s-%s" % (ch.strip() or "-", rs, rn, a, b), (rn, a, b)))
            if i + 1 < len(lst):
                nxt = lst[i + 1][1]["atoms"]
                if "C" in at and "N" in nxt and dist(at["C"], nxt["N"]) < 2.0:
                    out.append((at["C"], nxt["N"],
                                "%s%d-%d 肽键 C-N" % (ch.strip() or "-", rs, lst[i+1][0][1]),
                                ("*PEP*", "C", "N")))
    # 二硫键也是共价键 —— 不加进键表的话,check_clashes 会把每一对二硫键
    # 都报成"原子重叠"(SG-SG 约 2.03 Å,低于 2.2 Å 截断)。
    sgs = [((ck, rs), r["atoms"]["SG"]) for (ck, rs, ic), r in st.residues.items()
           if std(r["name"]) == "CYS" and "SG" in r["atoms"]]
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            if dist(sgs[i][1], sgs[j][1]) < 2.5:
                out.append((sgs[i][1], sgs[j][1],
                            "%s%d-%s%d 二硫键 SG-SG" % (sgs[i][0][0] or "-", sgs[i][0][1],
                                                     sgs[j][0][0] or "-", sgs[j][0][1]),
                            ("*SS*", "SG", "SG")))
    return out


def check_bond_lengths(st, rep):
    """B1 键长异常 —— 极小化能修,但会让建库工具报天文数字能量"""
    bad = []
    for (p, q, tag, (rn, a, b)) in _bond_list(st):
        d = dist(p, q)
        tol = BOND_TOL
        if rn == "*PEP*":
            ideal = PEPTIDE_CN
        elif rn == "*SS*":
            ideal, tol = 2.03, 0.20          # 二硫键
        else:
            ideal = PATCH_IDEAL.get(frozenset((a, b)))
            if ideal is None:
                ideal = IDEAL_SC.get((rn, frozenset((a, b))))
            if ideal is None:
                ideal = IDEAL_BONDS.get((a, b), IDEAL_BONDS.get((b, a)))
            if ideal is None:
                ideal = IDEAL_SIDECHAIN
        if abs(d - ideal) > tol:
            bad.append((abs(d - ideal), "%s = %.2f Å (应 %.2f)" % (tag, d, ideal)))
    if bad:
        bad.sort(reverse=True)
        rep.add("B", "FAIL", "键长异常 %d 处" % len(bad),
                "\n    ".join(t for _, t in bad[:12]) + ("\n    ..." if len(bad) > 12 else "") +
                "\n    能量极小化会修好这些,不影响最终结果 —— 但建库时单点能量会很高,属正常。")
    else:
        rep.add("B", "PASS", "键长正常")


def check_clashes(st, rep):
    """B2 原子重叠(网格近邻搜索,不依赖 numpy)"""
    pts, labels, owner = [], [], []
    for (ch, rs, ic), r in st.residues.items():
        # 非标准残基/配体/糖没有连接性模板,它们**内部**的共价键会被误当成重叠。
        # 但不能因此把整个残基从检测里剔除 —— 那样配体或封端基团与蛋白之间的
        # 真实冲突(实测可低至 0.4 Å)就完全看不见了。做法:保留原子,只排除
        # 这类残基的**残基内**原子对。
        known = std(r["name"]) in TEMPLATES
        for n, x in r["atoms"].items():
            pts.append(x)
            labels.append("%s%d%s:%s" % (ch.strip() or "-", rs, std(r["name"]), n))
            owner.append(None if known else (ch, rs, ic))
    bonded = set()
    idx = {id(p): i for i, p in enumerate(pts)}
    for (p, q, _t, _bt) in _bond_list(st):
        if id(p) in idx and id(q) in idx:
            bonded.add((min(idx[id(p)], idx[id(q)]), max(idx[id(p)], idx[id(q)])))
    # 1-3 邻居(共享一个成键原子)也排除
    nb = defaultdict(set)
    for i, j in bonded:
        nb[i].add(j); nb[j].add(i)
    excl = set(bonded)
    for i in nb:
        for j in nb[i]:
            for k in nb[j]:
                if k != i:
                    excl.add((min(i, k), max(i, k)))

    CUT = 2.2
    grid = defaultdict(list)
    for i, p in enumerate(pts):
        grid[(int(p[0]//CUT), int(p[1]//CUT), int(p[2]//CUT))].append(i)
    hits = []
    for cell, members in grid.items():
        neigh = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neigh.extend(grid.get((cell[0]+dx, cell[1]+dy, cell[2]+dz), []))
        for i in members:
            for j in neigh:
                if j <= i:
                    continue
                if (i, j) in excl:
                    continue
                # 非标准残基的残基内原子对:没有连接表,共价键会被误判成重叠
                if owner[i] is not None and owner[i] == owner[j]:
                    continue
                d = dist(pts[i], pts[j])
                if d < CUT:
                    hits.append((d, "%s <-> %s = %.2f Å" % (labels[i], labels[j], d)))
    if hits:
        hits.sort()
        rep.add("B", "FAIL", "原子重叠 %d 对 (< 2.2 Å 且非成键)" % len(hits),
                "\n    ".join(t for _, t in hits[:12]) + ("\n    ..." if len(hits) > 12 else "") +
                "\n    极小化通常能推开。若最短距离 < 1.0 Å,建库工具的单点能量会爆到千万级。")
    else:
        rep.add("B", "PASS", "无原子重叠")


def check_chain_breaks(st, rep):
    """B3 主链断裂 / 残基编号跳号"""
    breaks, jumps = [], []
    for ch, lst in st.chains().items():
        for i in range(len(lst) - 1):
            (_, rs0, _), r0 = lst[i]
            (_, rs1, _), r1 = lst[i + 1]
            if "CA" in r0["atoms"] and "CA" in r1["atoms"]:
                d = dist(r0["atoms"]["CA"], r1["atoms"]["CA"])
                if d > 4.2:
                    breaks.append("%s %d-%d (CA-CA %.2f Å)" % (ch.strip() or "-", rs0, rs1, d))
            if rs1 != rs0 + 1:
                jumps.append("%s %d -> %d" % (ch.strip() or "-", rs0, rs1))
    if breaks:
        rep.add("B", "FAIL", "主链断裂 %d 处" % len(breaks), ", ".join(breaks) +
                "\n    CA-CA 正常 3.8 Å。断裂说明中间缺残基,建库时会被当成两条独立链(各自加端基)。")
    else:
        rep.add("B", "PASS", "主链连续")
    if jumps:
        rep.add("C", "WARN", "残基编号不连续 %d 处" % len(jumps), ", ".join(jumps[:10]) +
                "\n    若结构上是连续的,只是编号跳号,通常无害;但要确认不是真的缺残基。")


def check_ramachandran(st, rep):
    """B4 主链二面角离群(宽松判据,只抓明显离谱的)"""
    out = []
    for ch, lst in st.chains().items():
        for i in range(1, len(lst) - 1):
            a0, a1, a2 = lst[i-1][1]["atoms"], lst[i][1]["atoms"], lst[i+1][1]["atoms"]
            rn = std(lst[i][1]["name"])
            if rn == "GLY":
                continue
            if not all(k in a1 for k in ("N", "CA", "C")) or "C" not in a0 or "N" not in a2:
                continue
            # 断链两侧算出的 φ/ψ 没有意义,先确认前后确实通过肽键相连
            if dist(a0["C"], a1["N"]) > 2.0 or dist(a1["C"], a2["N"]) > 2.0:
                continue
            phi = dihedral(a0["C"], a1["N"], a1["CA"], a1["C"])
            psi = dihedral(a1["N"], a1["CA"], a1["C"], a2["N"])
            ok = ((-180 <= phi <= -20 and (-90 <= psi <= 60 or 90 <= psi <= 180 or psi <= -150))
                  or (20 <= phi <= 100 and -20 <= psi <= 100))       # 左手螺旋区
            if not ok:
                out.append("%s%d%s (φ=%.0f ψ=%.0f)" % (ch.strip() or "-", lst[i][0][1], rn, phi, psi))
    if out:
        rep.add("B", "WARN", "主链二面角离群 %d 个" % len(out),
                ", ".join(out[:10]) + ("..." if len(out) > 10 else "") +
                "\n    落在 Ramachandran 禁区。平衡阶段通常会自己走回来,但数量多说明建模质量差。")
    else:
        rep.add("B", "PASS", "主链二面角正常")


# ── C 级 ─────────────────────────────────────────────────────────────────────

def check_integrity(st, rep):
    """解析完整性 —— 对 pdb 和 gro 都要跑。这是防"静默丢结构"的最后一道闸。"""
    if st.gro_split:
        rep.add("C", "WARN", "GRO 文件无链信息,已按 CA-CA 间距自动切分为 %d 条链" % st.gro_split,
                "若实际不是这么多条链(例如中间真的缺了残基),请改用带链信息的 PDB 复核。")
    if st.lost_atoms:
        rep.add("C", "FAIL", "有 %d 个原子在分组时被覆盖丢失" % st.lost_atoms,
                "说明多条链共用同一个链标识、且残基编号重复,后写入的把先写入的盖掉了。\n"
                "    PDB 请补 chain ID 或 segid;GRO 本身不带链信息,若自动切分失败请先拆成单链。\n"
                "    ⚠️ 不修的话,后面所有检查都只作用在最后一条链上,其余链等于没检查。")


def check_format(st, rep):
    if st.fmt != "pdb":
        return
    no_chain = [a for a in st.atoms if a.chain == " "]
    no_elem  = [a for a in st.atoms if not a.element]
    altloc   = [a for a in st.atoms if a.altloc not in (" ", "A")]
    occ      = [a for a in st.atoms if abs(a.occ - 1.0) > 1e-6 and a.occ != 0.0]
    icodes   = [a for a in st.atoms if a.icode != " "]
    short    = [(i, l) for i, l in st.raw_lines
                if l.startswith(("ATOM", "HETATM")) and len(l) < 78]

    if no_chain:
        segids = sorted(set(a.segid for a in no_chain if a.segid))
        rep.add("C", "FAIL", "第 22 列 chain ID 为空 (%d 个原子)" % len(no_chain),
                ("该文件是 CHARMM 风格:靠 73-76 列的 segid %s 标识链,但标准 PDB 解析器读第 22 列。\n"
                 "    CHARMM-GUI 等工具会直接报 'expected chain ID at column 22' 拒绝上传。\n"
                 "    修法:给每条链的第 22 列填一个字母。" % (segids or "(无)")))
    else:
        rep.add("C", "PASS", "chain ID 完整")

    if no_elem:
        rep.add("C", "WARN", "第 77-78 列元素符号缺失 (%d 个原子)" % len(no_elem),
                "多数工具能从原子名推断,但不是标准 PDB。建议补上。")
    if short:
        rep.add("C", "WARN", "%d 行短于 78 字符" % len(short),
                "行 %s ... —— 说明缺 segid/元素列。" % ", ".join(str(i) for i, _ in short[:5]))
    if altloc:
        # 分组时已经按残基选定了单一构象,所以这不再是阻塞性问题(否则几乎每个
        # 高分辨晶体结构都会被拦下),降为提示。
        rep.add("C", "WARN", "存在 altLoc 多重占位 (%d 个原子)" % len(altloc),
                "构象 %s。已自动只保留一套(优先主构象),请确认选的是你要的那一套。"
                % ", ".join(sorted(set(a.altloc for a in altloc))))
    if occ:
        rep.add("C", "WARN", "occupancy != 1 的原子 %d 个" % len(occ),
                "部分占位原子,通常来自晶体结构的无序区域,需要确认取舍。")
    if icodes:
        rep.add("C", "WARN", "存在插入码 (%d 个原子)" % len(icodes),
                "插入码 %s。部分工具处理不好,建议重新编号。"
                % ", ".join(sorted(set(a.icode for a in icodes))))

    # 同一残基里出现同名原子两次
    dup = []
    seen = defaultdict(set)
    for a in st.atoms:
        if a.altloc not in (" ", "", "A", "1"):
            continue                      # 次要构象本就重名,不算重复
        key = (chain_key(a), a.resseq, a.icode)
        if a.name in seen[key]:
            dup.append("%s%d:%s" % (chain_key(a) or "-", a.resseq, a.name))
        seen[key].add(a.name)
    if dup:
        rep.add("C", "FAIL", "重复原子 %d 个" % len(dup), ", ".join(dup[:10]))

    if st.extra_models:
        rep.add("C", "WARN", "文件含多个 MODEL(NMR ensemble),只检查了第一个模型",
                "若要检查其他模型,请先拆分文件。")
    if st.bad_coord_lines:
        rep.add("C", "FAIL", "坐标列无法解析的行 %d 行" % len(st.bad_coord_lines),
                "行 %s。多为坐标溢出(如 -1234.567 占满 9 列)或行被截断,这些原子已被丢弃。"
                % ", ".join(str(i) for i in st.bad_coord_lines[:8]))
    if st.dropped_altloc:
        rep.add("C", "WARN", "已丢弃次要构象原子 %d 个" % st.dropped_altloc,
                "只保留了主构象(altLoc 为空或 A)。请确认这是你想要的那一个。")

    has_end = any(l.startswith("END") for _, l in st.raw_lines)
    if not has_end:
        rep.add("C", "WARN", "缺 END 记录", "多数解析器容忍,但不规范。")

    if st.has_hydrogen:
        rep.add("C", "PASS", "文件含氢原子",
                "建库工具通常会删掉重建。若想保留,记得勾 'preserve hydrogen coordinates'。")


# ── D 级:膜蛋白 ──────────────────────────────────────────────────────────────

def check_membrane(st, rep, axis_idx, memb_half=20.0):
    """D 级 —— 需要 --membrane 开启。memb_half = 膜疏水半厚(Å),默认 20(P-P ~40 Å)"""
    ax_name = "xyz"[axis_idx]
    for ch, lst in st.chains().items():
        # names 与 cas 必须一一对应,否则后面按下标取会错位/越界
        pairs = [(std(r["name"]), r["atoms"]["CA"]) for _, r in lst if "CA" in r["atoms"]]
        if len(pairs) < 5:
            continue
        names = [p[0] for p in pairs]
        cas = [p[1] for p in pairs]
        tag = "链 %s" % (ch.strip() or "-")

        # D1 主轴与膜法向夹角
        v, c = principal_axis(cas)
        tilt = math.degrees(math.acos(min(1.0, abs(v[axis_idx]))))
        rep.add("D", "PASS", "%s 主轴与膜法向(%s)夹角 = %.1f°" % (tag, ax_name, tilt),
                "跨膜螺旋的天然倾角通常 10-30°。0° 说明是人为摆正的理想起点(平衡后会自己倾斜)。")

        # D2 拓扑朝向
        zN, zC = cas[0][axis_idx], cas[-1][axis_idx]
        rep.add("D", "PASS", "%s 拓扑:N 端 %s=%.1f,C 端 %s=%.1f -> N→C 指向 %s%s"
                % (tag, ax_name, zN, ax_name, zC, ax_name, "+" if zC > zN else "-"),
                "务必核对这与真实拓扑一致(如 I 型膜蛋白 C 端在胞质侧,对应含 PS/PE 的内叶)。")

        # D3 疏水段识别 + 疏水失配
        w, best, bi = 19, -1e9, 0
        if len(names) >= w:
            for i in range(len(names) - w + 1):
                s = sum(KD.get(n, 0.0) for n in names[i:i+w]) / w
                if s > best:
                    best, bi = s, i
            hyd_len = dist(cas[bi], cas[min(bi+w-1, len(cas)-1)])
            need = math.degrees(math.acos(min(1.0, 2*memb_half / hyd_len))) if hyd_len > 2*memb_half else 0.0
            rep.add("D", "PASS", "%s 最疏水 19 残基段 = 第 %d-%d 位 (KD 均值 %.2f),跨度 %.1f Å"
                    % (tag, bi+1, bi+w, best, hyd_len),
                    "膜疏水厚度按 %.0f Å 计:疏水失配需要约 %.0f° 倾角来补偿。"
                    % (2*memb_half, need))
            if best < 1.0:
                rep.add("D", "WARN", "%s 疏水性偏低 (最佳窗口 KD=%.2f)" % (tag, best),
                        "典型跨膜段 KD 均值 > 1.6。偏低说明这一段未必能稳定跨膜。")

        # D4 positive-inside:K/R 应集中在同一侧(胞质侧)
        pos = [(cas[i][axis_idx], names[i]) for i in range(len(names))
               if names[i] in ("LYS", "ARG")]
        if pos:
            mid = sum(p[axis_idx] for p in cas) / len(cas)
            up = sum(1 for z, _ in pos if z > mid)
            lo = len(pos) - up
            rep.add("D", "PASS", "%s 碱性残基 (Lys/Arg) 分布:%s 正向 %d 个,负向 %d 个"
                    % (tag, ax_name, up, lo),
                    "positive-inside 规则:K/R 应富集在胞质侧。确认这一侧放的是含 PS/PE 的内叶。")

        # D5 芳香带:Trp/Tyr 天然富集在界面区
        aro = [(cas[i][axis_idx], names[i]) for i in range(len(names)) if names[i] in ("TRP", "TYR")]
        if aro:
            mid = sum(p[axis_idx] for p in cas) / len(cas)
            at_if = sum(1 for z, _ in aro if abs(z - mid) > memb_half * 0.6)
            rep.add("D", "PASS", "%s 芳香残基 (Trp/Tyr) %d 个,其中 %d 个位于界面区"
                    % (tag, len(aro), at_if),
                    "天然跨膜蛋白的 Trp/Tyr 富集在膜-水界面。全在膜心说明取向可能不对。")


# ── E 级:多体系一致性 ────────────────────────────────────────────────────────

def check_consistency(structures, axis_idx):
    """E 级 —— 单看每个都正常、放一起才发现不对的那类问题"""
    if len(structures) < 2:
        return None
    rows, issues = [], []
    for st in structures:
        for ch, lst in st.chains().items():
            names = [std(r["name"]) for _, r in lst]
            cas = [r["atoms"]["CA"] for _, r in lst if "CA" in r["atoms"]]
            q = sum(res_charge(r["name"]) for _, r in lst)     # 用原始残基名,保留质子化态
            direction = "?"
            if len(cas) >= 2:
                direction = "+" if cas[-1][axis_idx] > cas[0][axis_idx] else "-"
            rows.append({
                "file": st.label, "chain": ch.strip() or "-", "n": len(names),
                "first": names[0] if names else "?", "last": names[-1] if names else "?",
                "charge": q, "dir": direction,
                "span": (max(p[axis_idx] for p in cas) - min(p[axis_idx] for p in cas)) if cas else 0.0,
            })
    # "?" 表示该链太短、无法判定朝向,不能当成"第三种朝向"参与比较,
    # 否则任何含短链/纯 HETATM 链的文件都会被误报成"朝向不一致"。
    dirs = set(r["dir"] for r in rows if r["dir"] != "?")
    if len(dirs) > 1:
        grp = defaultdict(list)
        for r in rows:
            if r["dir"] == "?":
                continue
            grp[r["dir"]].append("%s:%s" % (r["file"], r["chain"]))
        issues.append(("FAIL", "拓扑朝向不一致",
                       "  ".join("%s -> %s" % (k, ", ".join(v)) for k, v in grp.items()) +
                       "\n    这些体系要横向比较,朝向必须一致 —— 否则同一个残基面对的是不同叶层的脂质。"))
    lasts = set(r["last"] for r in rows)
    if len(lasts) > 1:
        issues.append(("WARN", "C 端残基不一致", "出现 %s。若构建时端基处理不同,能量不可比。" % ", ".join(sorted(lasts))))
    spans = [r["span"] for r in rows if r["span"] > 0]
    if spans and max(spans) - min(spans) > 8.0:
        issues.append(("WARN", "跨度差异较大",
                       "%.1f ~ %.1f Å。长度差异会影响疏水失配和倾角,比较时需要考虑。" % (min(spans), max(spans))))
    return rows, issues


# ─────────────────────────────────────────────────────────────────────────────
# 输出 / Reporting
# ─────────────────────────────────────────────────────────────────────────────

TIER_TITLE = OrderedDict([
    ("A", "A 级 — MD 永远修不好,建模前必须处理"),
    ("B", "B 级 — 极小化能修好(会让建库单点能量偏高,属正常)"),
    ("C", "C 级 — 文件格式,会导致上传/解析失败"),
    ("D", "D 级 — 膜蛋白拓扑与几何(需人工确认)"),
])
MARK = {"FAIL": "[✗]", "WARN": "[!]", "PASS": "[✓]"}


def print_report(rep, quiet):
    print("\n" + "=" * 78)
    print("  %s" % rep.name)
    print("=" * 78)
    for tier, title in TIER_TITLE.items():
        items = rep.items.get(tier)
        if not items:
            continue
        # D 级(膜蛋白)全部用 PASS 发出,但它们本质是"需人工确认的事实",不是通过项。
        # --quiet 下若一并隐藏,倾角/拓扑朝向/疏水失配/positive-inside 就全看不到了 ——
        # 而那恰恰是这个工具最需要人过目的部分。
        shown = [x for x in items if not (quiet and x[0] == "PASS" and tier != "D")]
        if not shown:
            continue
        print("\n%s" % title)
        print("-" * 78)
        for status, head, detail in shown:
            print("  %s %s" % (MARK[status], head))
            if detail and status != "PASS":
                for line in detail.split("\n"):
                    print("      %s" % line.strip() if not line.startswith("    ") else "  " + line)
            elif detail and not quiet:
                for line in detail.split("\n"):
                    print("      %s" % line.strip())


def main():
    ap = argparse.ArgumentParser(description="MD 建模前结构体检(零依赖)")
    ap.add_argument("files", nargs="+", help="PDB 或 GRO 文件")
    ap.add_argument("--membrane", action="store_true", help="打开 D 级膜蛋白检查")
    ap.add_argument("--normal", default="z", choices=["x", "y", "z"], help="膜法向轴,默认 z")
    ap.add_argument("--memb-half", type=float, default=20.0, help="膜疏水半厚 Å,默认 20")
    ap.add_argument("--quiet", action="store_true", help="只显示问题,隐藏通过项")
    args = ap.parse_args()

    axis_idx = {"x": 0, "y": 1, "z": 2}[args.normal]
    structures, worst = [], {"A": 0, "C": 0}

    for path in args.files:
        if not os.path.exists(path):
            print("跳过(文件不存在): %s" % path, file=sys.stderr)
            worst["C"] += 1
            continue
        try:
            st = Structure(path)
        except Exception as e:
            print("解析失败 %s: %s" % (path, e), file=sys.stderr)
            worst["C"] += 1
            continue
        if not st.residues:
            print("跳过(未找到蛋白质残基): %s" % path, file=sys.stderr)
            worst["C"] += 1
            continue
        structures.append(st)

        rep = Report("%s  —  %d 残基, %d 重原子, %d 条链"
                     % (st.name, len(st.residues),
                        sum(len(r["atoms"]) for r in st.residues.values()), len(st.chains())))
        check_chirality(st, rep)
        check_cis_peptide(st, rep)
        check_missing_atoms(st, rep)
        check_planarity(st, rep)
        check_ring_piercing(st, rep)
        check_disulfide(st, rep)
        check_bond_lengths(st, rep)
        check_clashes(st, rep)
        check_chain_breaks(st, rep)
        check_ramachandran(st, rep)
        check_integrity(st, rep)
        check_format(st, rep)
        if args.membrane:
            check_membrane(st, rep, axis_idx, args.memb_half)

        print_report(rep, args.quiet)
        worst["A"] += rep.counts("A")[0]
        worst["C"] += rep.counts("C")[0]

    # ── E 级 ──
    res = check_consistency(structures, axis_idx)
    if res:
        rows, issues = res
        print("\n" + "=" * 78)
        print("  E 级 — 多体系一致性对比")
        print("=" * 78)
        print("\n  %-28s %-6s %5s %6s %6s %7s %8s" %
              ("文件", "链", "残基", "首", "末", "净电荷", "跨度Å"))
        print("  " + "-" * 74)
        for r in rows:
            print("  %-28s %-6s %5d %6s %6s %+7d %8.1f" %
                  (r["file"][:28], r["chain"], r["n"], r["first"], r["last"], r["charge"], r["span"]))
        if issues:
            print()
            for status, head, detail in issues:
                print("  %s %s" % (MARK[status], head))
                for line in detail.split("\n"):
                    print("      %s" % line.strip())
                if status == "FAIL":
                    worst["A"] += 1
        else:
            print("\n  %s 各体系之间一致" % MARK["PASS"])

    # ── 总结 ──
    print("\n" + "=" * 78)
    if not structures:
        # "没东西可检"和"检过没问题"必须给不同的结论和退出码,
        # 否则接进流水线时会静默放行。
        print("  ✗ 没有成功解析出任何蛋白质结构 —— 这不等于检查通过。")
        print("=" * 78 + "\n")
        return 2
    if worst["A"] == 0 and worst["C"] == 0:
        print("  结论:未发现阻塞性问题,可以进入建模。")
    else:
        if worst["A"]:
            print("  ✗ 有 %d 项 A 级问题 —— 跑 MD 修不好,必须先处理。" % worst["A"])
        if worst["C"]:
            print("  ✗ 有 %d 项 C 级问题 —— 会导致上传/解析失败。" % worst["C"])
    print("=" * 78 + "\n")

    return (1 if worst["A"] else 0) | (2 if worst["C"] else 0)


if __name__ == "__main__":
    sys.exit(main())
