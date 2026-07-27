#!/usr/bin/env python3
"""CG-MD (MARTINI) 自发组装体系的建库前体检 — 零依赖,纯 python3 标准库。

在把一个 MARTINI 蛋白 .itp / production .mdp 送去跑几十微秒之前,检查那几个"一旦错了
就白烧 GPU-天"的致命点。这些点全都能从纯文本文件里静态判定,不需要 GROMACS/numpy。

用法:
    cg_preflight.py PROA_0.itp [more.itp ...]          # 检查拓扑
    cg_preflight.py prod.mdp                            # 检查 mdp
    cg_preflight.py PROA_0.itp prod.mdp                # 混合,自动按扩展名分派

退出码:0 = 无阻塞问题;1 = 有 FATAL(跑之前必须修)。
"""
import sys
import re

# ---- 分级 --------------------------------------------------------------------
# FATAL : 跑之前必须修,否则要么秒崩,要么结果从物理上就是错的
# WARN  : 大概率是问题,但取决于意图,人工确认
# INFO  : 只报事实(净电荷等),供你核对
FATAL, WARN, INFO = "FATAL", "WARN", "INFO"
ICON = {FATAL: "✗ FATAL", WARN: "! WARN ", INFO: "· info "}


def _sections(lines):
    """把 .itp/.top 切成 {section_name: [(lineno, raw), ...]},去注释、去空行。"""
    secs, cur = {}, None
    for i, raw in enumerate(lines, 1):
        line = raw.split(";", 1)[0].rstrip()
        if not line.strip():
            continue
        m = re.match(r"\s*\[\s*(\S+)\s*\]", line)
        if m:
            cur = m.group(1).lower()
            secs.setdefault(cur, [])
            continue
        if cur is not None:
            secs[cur].append((i, line))
    return secs


def check_itp(path):
    """检查单个 MARTINI 蛋白 moleculetype .itp。返回 (findings, netcharge)."""
    out = []
    with open(path) as fh:
        lines = fh.readlines()
    secs = _sections(lines)

    # --- 1. 刚性键软化(fix_stiff_bonds 必做) ---------------------------------
    # martinize2 把环/侧链刚性连接写成 k=1e6 的 [bonds](振动周期~30 fs)→ dt≥15 fs 必炸。
    # fix_stiff_bonds.py 应已把它们软化到 k=1e5。这里查 [bonds] 里还有没有漏网的 k≥5e5。
    stiff = []
    for ln, row in secs.get("bonds", []):
        c = row.split()
        if len(c) >= 5:
            try:
                k = float(c[4])
            except ValueError:
                continue
            if k >= 5e5:
                stiff.append((ln, k))
    if stiff:
        out.append((FATAL,
            f"[bonds] 里有 {len(stiff)} 条 k≥5e5 的刚性键(如 line {stiff[0][0]}: k={stiff[0][1]:.0f})"
            " —— fix_stiff_bonds.py 没跑或没生效。dt≥15fs(甚至 8fs)会炸。"
            " 修:对该 itp 跑 01_martinize/fix_stiff_bonds.py 软化到 k=1e5。"))
    elif "bonds" in secs:
        out.append((INFO, "[bonds] 无 k≥5e5 刚性键 —— fix_stiff_bonds 已生效。"))

    # --- 2. 弹性网络存在,且只在链内(martinize per-monomer 的前提) ------------
    # 这个 itp 是单个 moleculetype;只要弹性网络写在这个 moleculetype 内、原子序号都在本链范围,
    # 就是链内的。真正致命的"跨链弹性网络"发生在把整个三聚体一次性 martinize 时。
    # 这里能查的是:(a) 有没有弹性网络;(b) itp 里只有一个 [moleculetype](不是多链拼在一起)。
    n_moltype = len(secs.get("moleculetype", []))
    if n_moltype > 1:
        out.append((FATAL,
            f"这个 itp 含 {n_moltype} 个 [moleculetype] —— 像是多条链拼在一个 itp 里 martinize 的产物,"
            " 弹性网络/exclusions 很可能跨了链,会把三聚体人为焊死。"
            " 修:对**单个 monomer** 分别 martinize,每链一个 itp。"))
    has_en = any("rubber" in r.lower() or "elastic" in r.lower()
                 for r in [l for _, l in lines_with_comments(path)])
    # rubber band 注释被 _sections 去掉了,单独扫一遍原始行找标记 + 找长程 bond
    if not has_en:
        # 退而求其次:弹性网络表现为跨多个残基的 [bonds](i,j 相隔较远、k~500)
        longrange = 0
        for _, row in secs.get("bonds", []):
            c = row.split()
            if len(c) >= 5:
                try:
                    i, j, k = int(c[0]), int(c[1]), float(c[4])
                except ValueError:
                    continue
                if abs(i - j) >= 6 and 100 <= k <= 2000:
                    longrange += 1
        has_en = longrange >= 3
    if has_en:
        out.append((INFO, "检出链内弹性网络(elastic network)—— 单链结构会被维持。"))
    else:
        out.append((WARN,
            "没检出弹性网络。TM 短螺旋没有 EN 在长模拟里可能解折叠。"
            " 确认 martinize2 用了 -elastic(本项目 -ef 500 -el 0.5 -eu 0.9)。"))

    # --- 3. 净电荷(中性帽 + 体系中和的核对值) --------------------------------
    q = 0.0
    n_atoms = 0
    for _, row in secs.get("atoms", []):
        c = row.split()
        if len(c) >= 7:
            try:
                q += float(c[6]); n_atoms += 1
            except ValueError:
                pass
    qr = round(q)
    out.append((INFO,
        f"单链净电荷 ≈ {qr:+d} e({n_atoms} 珠)。三拷贝体系需 genion 中和到 0.15 M NaCl 时用这个真值"
        "(insane 读不出 CG 蛋白电荷)。若 |小数部分| 偏离整数较多,检查帽子/质子化态。"))
    # 中性帽:N/C 端 BB 珠不应带 ±1(那是没加中性帽的短肽假电荷)
    atoms = secs.get("atoms", [])
    if atoms:
        first_bb = next((r for _, r in atoms if len(r.split()) >= 5 and r.split()[4] == "BB"), None)
        last_bb = next((r for _, r in reversed(atoms) if len(r.split()) >= 5 and r.split()[4] == "BB"), None)
        for tag, r in (("N 端", first_bb), ("C 端", last_bb)):
            if r:
                cc = r.split()
                try:
                    if abs(float(cc[6])) >= 0.5:
                        out.append((WARN,
                            f"{tag} BB 珠带电 {float(cc[6]):+.1f} —— 可能没加中性帽(-nt)。"
                            " 孤立短肽的末端 ±1 电荷在真实全长蛋白里不存在,会污染组装。"))
                except (IndexError, ValueError):
                    pass
    return out, qr


def lines_with_comments(path):
    with open(path) as fh:
        return list(enumerate(fh.readlines(), 1))


def check_mdp(path):
    """检查 MARTINI production/equil .mdp。"""
    out = []
    kv = {}
    for raw in open(path):
        line = raw.split(";", 1)[0]
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip().lower().replace("_", "-")] = v.strip().lower()

    dt = kv.get("dt")
    if dt is not None:
        try:
            if float(dt) > 0.008 + 1e-9:
                out.append((FATAL,
                    f"dt = {dt} > 0.008 —— MARTINI 本项目实测:20fs ~20ps 崩、15fs ~150ps 崩、"
                    "10fs 长生产偶发 LJ 硬碰撞崩。生产用 8 fs(dt=0.008)。"))
            else:
                out.append((INFO, f"dt = {dt}(≤0.008,OK)。"))
        except ValueError:
            pass

    ct = kv.get("coulombtype", "")
    if "pme" in ct:
        out.append((FATAL,
            f"coulombtype = {ct} —— MARTINI 用 reaction-field,不是 PME。"
            " 且 mdrun 绝不能加 -pme gpu(fatal)。改 coulombtype = reaction-field。"))
    elif "reaction-field" in ct or "reaction_field" in ct:
        out.append((INFO, "coulombtype = reaction-field(OK,记得 mdrun 只 -nb gpu、不要 -pme gpu)。"))
    elif ct:
        out.append((WARN, f"coulombtype = {ct} —— MARTINI 通常用 reaction-field,确认这是有意的。"))

    er = kv.get("epsilon-r")
    if er and er not in ("15", "15.0"):
        out.append((WARN, f"epsilon-r = {er} —— MARTINI 标准是 15,确认。"))

    for key, want, note in [
        ("rvdw", "1.1", "MARTINI 3 标准截断"),
        ("rcoulomb", "1.1", "MARTINI 3 标准截断"),
    ]:
        v = kv.get(key)
        if v and v.rstrip("0").rstrip(".") != want.rstrip("0").rstrip("."):
            out.append((WARN, f"{key} = {v}(MARTINI 3 常用 {want};{note})。"))

    if not out:
        out.append((WARN, f"{path}:没识别到 dt/coulombtype 等关键项,确认这是个 mdp。"))
    return out


def main(argv):
    files = [a for a in argv[1:] if not a.startswith("-")]
    if not files:
        print(__doc__)
        return 0
    worst_fatal = False
    charges = {}
    for path in files:
        print(f"\n=== {path} ===")
        try:
            if path.endswith((".mdp",)):
                findings = check_mdp(path)
            elif path.endswith((".itp", ".top")):
                findings, q = check_itp(path)
                charges[path] = q
            else:
                # 按内容猜:含 [ moleculetype ] 当 itp,否则当 mdp
                head = open(path).read(4000).lower()
                if "[ moleculetype" in head or "[moleculetype" in head:
                    findings, q = check_itp(path); charges[path] = q
                else:
                    findings = check_mdp(path)
        except FileNotFoundError:
            print(f"  {ICON[FATAL]}  文件不存在"); worst_fatal = True; continue
        for level, msg in sorted(findings, key=lambda x: [FATAL, WARN, INFO].index(x[0])):
            print(f"  {ICON[level]}  {msg}")
            if level == FATAL:
                worst_fatal = True

    if charges:
        print("\n=== 净电荷汇总(genion 中和用真值)===")
        for p, q in charges.items():
            print(f"  {q:+d} e  x3 拷贝 = {3*q:+d} e 需中和   {p}")

    print()
    if worst_fatal:
        print("结论:有 FATAL —— 跑之前必须修,否则白烧 GPU-天。")
        return 1
    print("结论:无阻塞问题。WARN 项请对照意图人工确认。")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
