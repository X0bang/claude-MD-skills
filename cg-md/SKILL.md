---
name: cg-md
description: 用 MARTINI 3 粗粒化(CG)自发组装 MD 研究跨膜(TM)螺旋能否/多快自组装成寡聚体的完整工作流与避坑指南。当用户要做 CG-MD、MARTINI、martinize2/insane 建库、TM 螺旋自发组装/寡聚化、蛋白-蛋白侧向组装动力学,或要解读 CG 组装结果(first-passage time、组装概率、界面验伪)时使用。也用于判断 CG 结果可不可信、要不要用原子级方法交叉检验。
allowed-tools: Bash(${CLAUDE_SKILL_DIR}/cg_preflight.py *) Read Bash(module load *) Bash(gmx *) Bash(sbatch *) Bash(squeue *)
---

# MARTINI CG 自发组装 MD

用 MARTINI 3 把"几条已各自插膜的 TM 螺旋能不能自发侧向组装成寡聚体、多快、多大概率"直接跑出来。
**这个 skill 的核心价值有两条:(1) 那几个一旦错了就白烧 GPU-天的建库/跑法坑;(2) 怎么判断 CG 结果
可不可信 —— CG 组装强度排序经常是力场假象,必须用原子级方法(US/PMF)交叉检验。**

跑之前先体检:

```bash
python3 ${CLAUDE_SKILL_DIR}/cg_preflight.py PROA_0.itp prod.mdp
```

零依赖(纯 python3 标准库),静态查那几个致命点:刚性键有没有软化、弹性网络在不在/是不是跨链、
净电荷(genion 中和用)、末端中性帽、dt≤8fs、coulombtype=reaction-field。FATAL 项跑之前必须修。

---

## 五个绝不能违反的科学约束(错了结果从物理上就是错的)

这些是本方法成立的前提,违反任何一条,"组装与否"就不再由真实物理力决定:

1. **弹性网络(EN)只能加在单条链内部,绝不能跨链。** 做法上靠**对单个 monomer 分别 martinize**
   从构造上保证 EN 只在链内 —— 千万别把整个三聚体一次性 martinize(EN/exclusions 会跨链,
   等于把寡聚体人为焊死,组装就成了预设结论)。`cg_preflight.py` 会警告一个 itp 里出现多个
   `[moleculetype]`(多链拼装的迹象)。
2. **末端加中性帽(乙酰化/酰胺化,`martinize2 -nt`)。** 否则孤立短肽的 N/C 端带上真实全长蛋白里
   不存在的 ±1 电荷,污染组装能量学。
3. **二级结构必须手动喂给 martinize2。** 设计的短肽在 PDB 数据库里没有条目可推断螺旋 ——
   对 monomer 跑 DSSP 拿到 SS 字符串,`martinize2 -ss <字符串>`。
4. **膜配方对齐原子级建模。** 脂质种类/比例照搬 CHARMM-GUI 的 `[ molecules ]`(本项目
   POPC:POPE:POPS≈45:38:16),对称双层,0.15 M NaCl。否则膜环境不可比。
5. **远间距起始 + 周期镜像自检。** 三拷贝摆得足够远(质心距 > 盒边/2,使真实距=最小镜像距),
   建完**强制断言**两两最小镜像距 ≥ 阈值(本项目 7.5 nm),否则组装可能是镜像假象。

---

## 工作流(五步)

按**步骤**分目录,每步自包含(脚本 + 四变体产物 + 日志都在里面);变体作子目录横切,重复只在生产步出现。

1. **martinize** —— `martinize2`(`vermouth` 包)把每个 monomer.pdb 转 MARTINI 3。
   **逐 monomer**、`-elastic -ef 500 -el 0.5 -eu 0.9`、`-nt`(中性帽)、`-ss <DSSP 串>`。
   **然后必须跑 `fix_stiff_bonds`(见下,不可跳)。**
2. **build** —— PCA 把螺旋轴对齐 z → 三拷贝摆等边三角形(centroid 居中)→ `insane.py` 建膜
   (`-center` 让螺旋跨双层)→ `gmx genion` 用**真拓扑电荷**中和到 0.15 M NaCl
   (insane 读不出 CG 蛋白电荷,得自己算,`cg_preflight.py` 会报净电荷)→ 镜像距自检。
   **run 目录深度要固定**(topol.top 用相对路径 include 力场,本项目要求 3 层深)。
3. **equilibration** —— EM → 平衡阶梯(5fs→10fs,主链位置约束,末段释放)→ 出 `eq3.{gro,cpt}`。
   tcoupl 分三组 SOLU/MEMBRANE/SOLVENT(`gmx select` 建)。
4. **production** —— **dt=8 fs**、reaction-field(**不能 -pme gpu**)、`-nb gpu -update cpu`;
   多重复(靠初始方位角独立化)× 几十 µs;SLURM,可抢占 QOS 要能自我续投 + checkpoint 恢复。
5. **analysis** —— 逐帧**链间 BB 骨架珠最小镜像距离**(box-aware,无需 unwrap)→ **阈值 + 持续时间**
   判据(本项目 0.65 nm + ≥200 ns)→ 连通分量分级 dimer→trimer + first-passage time →
   组装概率;**界面验伪**:自组装界面 vs 设计静态界面(recovery/precision/jaccard),
   剔除 MARTINI 非特异过黏假阳性。

---

## 硬核坑(全都踩过,违反就重新崩)

| 坑 | 后果 | 解 |
|---|---|---|
| **不软化刚性键** | martinize2 把环/侧链刚性连接写成 k=1e6 *bond*(周期~30fs),dt≥15fs 必炸 | `fix_stiff_bonds.py` 软化到 **k=1e5**(周期~90fs,涨落~0.005nm,实质仍刚性)。**别转成 constraints**(过耦合主链约束会崩) |
| **dt 太大** | 20fs ~20ps 崩、15fs ~150ps 崩、10fs 长生产偶发 LJ 硬碰撞崩 | 生产 **8 fs**,把最爱崩的构型跑干净,只多花 ~20% 机时 |
| **`-pme gpu`** | MARTINI 用 reaction-field,加 PME GPU 直接 fatal | 只 `-nb gpu`,不要 `-pme gpu` |
| **GPU update** | 这套约束拓扑不能用 GPU update groups,GPU-resident update+constrain 在 step 0 死锁 | **`-update cpu`**(nonbonded 仍在 GPU) |
| **一卡塞多条** | 小体系 CPU-update 受限、latency-bound,聚合吞吐不升反降;卡型也几乎无差别(~2000–2870 ns/day) | **一卡一条**,靠多卡并行 |
| **可抢占 QOS 直接跑** | 被抢占就丢进度 | `--requeue` + `-cpt 3`(3分钟 checkpoint)+ `-cpi` 续跑 + 自我续投(把卡型/QOS 粘住带走) |
| **DSSP 版本** | vermouth `-dssp` 只认 v2/v3 CLI;v4 又拒绝空 chain-ID 的 protomer PDB | 绕:Biopython 盖 chain A → mmCIF → `mkdssp` → 把 SS 串喂 `martinize2 -ss`(见 `01_martinize/dssp_ss.py`) |
| **insane 报 pkg_resources** | insane 依赖已被移除的 `pkg_resources` | 环境钉 **setuptools<81**,别升级 |
| **sbatch spool 路径** | sbatch 把脚本拷到 spool,`BASH_SOURCE` 指向 spool,找不到项目根 | 用 `${SLURM_SUBMIT_DIR:-<BASH_SOURCE 兜底>}` 定位根 |
| **连续 resid** | 三拷贝在 tpr 里 resid 连续(1–30/31–60/61–90) | 逐残基界面占据率要把珠映射到**每链本地**残基位,别用原始 resid |

---

## 怎么读结果 —— 最重要的部分

### 组装判据本身
- **first-passage time / 组装概率**:稳健,可报。
- **`final_state`(末帧分类)脆弱**:链到达三聚体后可能又散开,末帧口径会**低估**组装。
  必须与 **"曾到达过(ever_trimer)"** 口径并列。(实例:某变体 3 条 replica 末帧只 1 条是三聚体,
  但 3 条**都**记录到了 trimer first-passage —— 实际是 3/3 曾组装。)

### 界面验伪(判特异 vs 过黏)
MARTINI 蛋白-蛋白相互作用偏"黏",会产生**非特异聚集假阳性**。所以**光看"组装没组装"不够**,
必看自组装界面是否落在设计界面上:
- **precision 高(如 0.93)** = 组装出来的接触面确实在设计界面上 → **特异组装,可信正面结论**。
- precision 低 = 非特异过黏,别当组装信号。

### ⚠️ CG 的跨变体"组装强度排序"经常不可信 —— 必须交叉检验
这是最容易被误用的一点。**MARTINI 4 重原子/珠的粗粒化会:**
- **欠稳定**特异堆积(大侧链 Leu/Ile/Val 的 knobs-into-holes、芳香 Phe/Tyr/Trp 堆叠);
- **相对高估**小残基/泛疏水的泛化接触(Ala 密集的平界面靠泛化 LJ 就搭得上)。

后果:一个界面"少而强"(芳香/大侧链、井深)的变体,CG 可能判它**不组装**;一个"多而浅"
(Ala 多、footprint 大但单点弱)的变体,CG 可能判它**强组装**——**顺序可以完全反**。

**判据:接触残基多 ≠ 结合强。footprint 大小不是结合自由能。** 结合强度以**原子级 PMF/自由能**为准。
**当 CG 排序与原子级 PMF 矛盾时,信原子级** —— CG 的贡献退回到"能不能观测到特异组装"这个定性、
单体系的结论,不是跨变体强度排序。

### 和湿实验对接前先分清两条正交轴
- **表达量 / 表面展示(如 FACS)← 插入 ΔG(逐螺旋、translocon topogenesis 门控)。**
  单体插入但不组装照样有信号;线性肽表位的抗体检测**对是否组装是盲的**。
- **组装 ← 界面互补性(另一条轴)。**
两轴无因果关系。用 CG 组装去"解释表达量"多半是牵强的事后合理化 —— 表达讲 ΔG,组装讲界面/PMF。
测三聚体本身要 FRET/BRET、BN-PAGE、交联、分子亮度法,不是普通 FACS。

---

## 用完之后
1. **改过 itp/mdp 一定重跑 `cg_preflight.py`**,确认 FATAL 清零(尤其别在改动里重新引入 k≥5e5 刚性键)。
2. **组装强度结论别只靠 CG** —— 若手上有原子级 US/PMF,拿来交叉检验;矛盾时以原子级为准,
   并在报告里把 CG 的偏差当方法学结论写清楚。
3. `final_state` 与 `ever_trimer` 两个口径并列报,别只报末帧。
