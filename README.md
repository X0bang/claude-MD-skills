# claude-MD-skills

分子动力学(MD)工作流用的 Claude Code skills。

## 已有 skills

| Skill | 作用 |
|---|---|
| [`md-structure-check`](md-structure-check/) | 跑 MD 之前对蛋白质结构做体检,按"MD 修得好 / 修不好"分级报告。零依赖,纯 python3 标准库。 |
| [`cg-md`](cg-md/) | MARTINI 3 粗粒化自发组装 MD 的完整工作流 + 避坑 + 结果解读(含 CG 强度排序不可信、须与原子级 PMF 交叉检验)。附零依赖建库前体检脚本 `cg_preflight.py`。 |

## 安装

Claude Code 从 `~/.claude/skills/<名字>/SKILL.md` 读取个人 skill,并且**支持符号链接**。
所以把仓库克隆到任意位置,再把需要的 skill 链进去即可:

```bash
git clone git@github.com:X0bang/claude-MD-skills.git ~/claude-MD-skills
mkdir -p ~/.claude/skills
ln -s ~/claude-MD-skills/md-structure-check ~/.claude/skills/md-structure-check
```

之后在任意项目里输入 `/md-structure-check` 调用,或者直接描述需求让 Claude 自动识别。

更新所有机器:

```bash
cd ~/claude-MD-skills && git pull
```

> ⚠️ 如果 `~/.claude/skills/` 这个目录在 Claude Code 启动时还不存在,新建之后需要**重启 Claude Code**
> 才能被识别。目录已存在的话,增删 skill 会实时生效。

## 脚本也可以单独用

不装 skill 也能当普通命令行工具跑:

```bash
python3 ~/claude-MD-skills/md-structure-check/check_structure.py protein.pdb --membrane
```

退出码:`0` 无阻塞问题 / `1` 有 A 级问题 / `2` 有 C 级问题 / `3` 两者都有 —— 可用于自动化流水线。

## 设计原则

**零依赖。** 集群、登录节点、别人的机器上未必有 numpy 或 conda 环境。所有脚本只用 python3 标准库。

**分级而非罗列。** 结构检查工具的常见毛病是把几十条问题平铺出来,让人分不清哪条要命。
这里强制区分"跑 MD 能自愈的"和"永远修不好的" —— 后者才是必须立刻处理的。
