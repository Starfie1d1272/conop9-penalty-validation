# CONOP 代价函数的 Python 复现与校验

**作者**：课程小组成员二  
**单位**：南京大学地球科学与工程学院，江苏 南京 210023

**摘要**：原版 CONOP9 程序能够完成高维地层序列优化，但其 Windows 二进制运行方式和封闭实现限制了自动化实验、误差诊断和算法扩展。本文基于 `geology-big-data-hw-main` 项目的 Python 复现代码，系统说明 CONOP 输入解析、复合序列表示、Ordinal、Level 与 Eventual 三类代价函数的实现逻辑，并以原版 `bestsoln.dat` 和 `outmain.txt` 为基准进行校验。项目中 `conop_py/io.py` 负责解析 `sections.txt`、`events.txt`、`loadfile.dat` 与 `bestsoln.dat`，`conop_py/cost.py` 通过 `ConopContext` 预计算剖面观测、事件位置和 taxon 延限信息。结果显示，Python 版 Ordinal misfit 在原版解上得到 367，与 CONOP9 完全一致；Level misfit 经 L1 保序回归和逐剖面校正后得到 237，与 CONOP9 输出一致；Eventual misfit 得到 335，较 CONOP9 的 353 低 18，差异主要可能来自 PAV 保序回归中偶数块合并的中位数选择细节。该复现工作不仅验证了项目核心数据流，也为后续快速模拟退火、多目标优化与不确定性分析提供了可解释、可测试的算法基础。

**关键词**：CONOP；Python 复现；代价函数；保序回归；Ordinal；Level；Eventual

**中图分类号**：P539；TP311.1  
**文献标识码**：A  
**文章编号**：课程论文-2026-02

## Python Reimplementation and Validation of CONOP Misfit Functions

**Abstract**: This paper describes a Python reimplementation of CONOP misfit functions based on the `geology-big-data-hw-main` project. The implementation parses CONOP input files, represents composite event sequences, and evaluates Ordinal, Level and Eventual penalties through a shared `ConopContext`. Validation against the original CONOP9 solution shows exact agreement for Ordinal and Level penalties, while Eventual differs by 18 points, probably due to tie-breaking details in isotonic regression. The reimplementation provides a transparent foundation for automated experiments and algorithmic extensions.

**Key words**: CONOP; Python reimplementation; misfit function; isotonic regression; Ordinal; Level; Eventual

## 1 引言

CONOP 的核心思想是将多剖面地层对比问题转化为事件排列问题：所有化石 FAD、LAD 和 marker 被放入一个全局复合序列中，算法寻找使各剖面观测冲突最小的排列。原版 CONOP9 程序在地层学研究中应用广泛，但其二进制程序不便于逐步调试、批量实验和新目标函数开发。对于课程项目而言，仅使用原版程序可以得到结果，却难以回答“为什么该解较好”“某一事件为什么不稳定”“是否能替换接受准则”等更深入的问题。

本项目因此开发了 Python 复现版。该复现不是简单调用 CONOP 输出，而是从输入数据解析、事件序列构建、代价函数计算到模拟退火搜索均重新实现。本文聚焦代价函数部分，说明项目如何将地层约束转化为可计算的 misfit，并讨论与 CONOP9 对齐过程中的关键问题。

## 2 数据结构与输入解析

### 2.1 输入文件解析

项目中 `conop_py/io.py` 定义了四类核心数据结构：`Section`、`Entity`、`Observation` 与 `SolutionRecord`。其中 `Observation` 保存剖面编号、实体编号、事件类型、层位以及权重信息，是构建代价函数的基础。

事件类型采用项目中确认的编码：

| 类型码 | 含义 | 说明 |
|---:|---|---|
| 1 | FAD | 分类单元首现 |
| 2 | LAD | 分类单元末现 |
| 4 | ASH | 火山灰层等 marker |
| 5 | AGE | 同位素年龄约束 marker |

解析过程特别注意 taxon 与 marker 的区分。`events.txt` 中实体编号并不总是简单连续，因此代码通过 `infer_taxa_from_observations()` 从 `loadfile.dat` 中出现 FAD 或 LAD 的实体推断 taxon 集合，再用 `parse_events()` 标记实体类型。这避免了按编号硬切分导致的错配。

### 2.2 复合序列表示

在 Python 版中，一个事件由二元组 `(entity_id, event_type)` 唯一表示。例如某菊石的 FAD 表示为 `(eid, 1)`，LAD 表示为 `(eid, 2)`。`solution_to_sequence()` 将 `bestsoln.dat` 中的三列记录按 position 排序，转换为事件二元组列表。所有代价函数都围绕这个列表计算事件在复合序列中的位置。

### 2.3 ConopContext 预计算结构

`conop_py/cost.py` 中的 `ConopContext` 是复现版的中心数据结构。它在给定模型序列和剖面观测后预先计算：

| 字段 | 作用 |
|---|---|
| `pos` | 事件到复合序列位置的映射 |
| `section_obs` | 每个剖面的观测事件列表 |
| `sec_levels` | 每个剖面按层位聚合的事件 |
| `taxon_sec` | 每个 taxon 在每个剖面中的 FAD/LAD 层位 |
| `taxa` | 同时具有 FAD 和 LAD 的分类单元集合 |

这种设计的意义在于，各类 misfit 共享同一套预计算信息。序列变化时只需用 `rebuild_pos()` 更新事件位置，而不必反复解析数据或重建剖面结构。后续模拟退火中的增量计算也建立在这一结构之上。

## 3 三类代价函数的实现

### 3.1 Ordinal misfit

Ordinal 惩罚统计剖面内观测顺序与复合序列顺序不一致的事件对数。实现时，代码对每个剖面按照实际层位排序，再将事件替换为其在复合序列中的 rank，最后用逆序对计数得到该剖面的 ordinal 罚分。所有剖面罚分相加得到总值。

在原版 `CONOP-run/bestsoln.dat` 上，Python 版 `ordinal_misfit()` 得到 367，与 CONOP9 输出完全一致。这说明输入解析、事件键表示和基本顺序约束均正确，是复现工作的第一道关键校验。

### 3.2 Level misfit

Level 惩罚的目标不是简单计算逆序对，而是衡量为了使某一 taxon 的延限范围与复合序列一致，需要在剖面内跨过多少个 distinct horizon。项目采用 L1 保序回归（PAV, pool adjacent violators）计算事件在剖面中的“放置水平”，再比较放置水平与观测水平之间跨越的 horizon 数。

代码中 `_compute_placed_isotonic()` 对剖面内事件按复合序列位置排序，并为不同类型事件设置 box constraints：

| 事件类型 | 放置约束 |
|---|---|
| FAD | 放置水平不高于观测 FAD |
| LAD | 放置水平不低于观测 LAD |
| marker | 固定在观测层位 |

若相邻块违反单调性，则 PAV 算法合并块，并用下中位数作为 L1 目标下的代表值。最终 `_sec_level()` 统计每个 FAD 或 LAD 需要向上或向下跨越的 horizon 数。逐剖面结果与 CONOP9 对齐后，`level_misfit()` 在原版解上得到 237，与 CONOP9 的 Level Penalty 一致。

### 3.3 Eventual misfit

Eventual 与 Level 的地质含义相近，但跨越某一 horizon 时不再只计 1 分，而是按该 horizon 上位于 taxon 复合延限内部的 forcing event 数加权。因此 Eventual 更强调被延限扩展影响的事件数量，较适合识别局部层位上事件聚集造成的冲突。

项目中 `_sec_eventual()` 复用保序回归放置水平，随后对跨越的每个 horizon 查找该层位上的事件数，并排除 AGE 类型事件。该实现得到 Python Eventual = 335，而 CONOP9 输出为 353，相差 18。由于 Level 已完全对齐，差异不来自输入解析或基本延限逻辑，而更可能来自 Eventual 计数细节或保序回归合并块的 tie-breaking。

## 4 校验结果与误差诊断

### 4.1 与 CONOP9 的总体对比

在同一个原版 `bestsoln.dat` 上，三类罚分对比如表 1。

**表 1 Python 复现版与 CONOP9 的代价函数对比**

| 指标 | Python | CONOP9 | 差异 | 评价 |
|---|---:|---:|---:|---|
| Ordinal | 367 | 367 | 0 | 完全一致 |
| Level | 237 | 237 | 0 | 完全一致 |
| Eventual | 335 | 353 | -18 | 存在小幅差异 |

Ordinal 和 Level 的完全一致，表明项目已经正确复现了 CONOP 最重要的顺序约束与 level 罚分。Eventual 的 -18 差异约占 CONOP9 Eventual 的 5.1%，在地层解释上不改变主要趋势，但对算法对齐具有诊断价值。

### 4.2 Eventual 差异的集中位置

项目中的 `scripts/eventual_diagnosis.py` 对 Eventual 贡献进行了剖面和事件拆分。按剖面看，贡献最大的剖面为 Seymour Island C、Seymour Island D、Seymour Island F、Seymour Island A 与 Quiriquina Island；按事件看，Top 10 事件累计贡献约 60%，达到 80% 累计贡献需要前 17 个事件。高贡献事件包括 `Kitchinites darwini` LAD、`Grossouvrites gemmatus` LAD、`Maorites seymourianus` LAD、`Zelandites varuna` FAD 等。

这说明 Eventual 差异不是均匀分散在所有事件上，而是集中于少数剖面和少数高冲突事件。对课程报告而言，这些事件可作为“矛盾热点”的实例，说明自动化代价函数不仅能给出总分，还能定位问题来源。

### 4.3 PAV 中位数选择的影响

项目报告提出的最可能解释是 PAV 保序回归在偶数块合并时的中位数选择差异。Python 版 `_lower_median()` 在偶数个观测值中取下中位数；而 CONOP9 的内部实现可能取上中位数或采用略不同的块合并顺序。若 5 至 15 个事件的放置水平因此相差 1 个 horizon，每个事件贡献 1 至 2 分，就可解释约 18 分差距。

这类差异具有方法论意义：对于封闭二进制程序，完全复现不只依赖数学定义，还依赖若干未公开的实现细节。Python 复现版的优势在于所有假设都显式写在代码中，可通过测试锁定，也可在后续实验中替换不同规则。

## 5 工程实现的进一步价值

### 5.1 支持模拟退火自动化

在代价函数可复现后，`conop_py/anneal.py` 实现了模拟退火主循环。`AnnealConfig` 包含 `startemp`、`ratio`、`steps`、`trials`、`seed`、`force_fad_before_lad`、`coex_penalty`、`accept_rule` 等配置项。相较原版 CONOP9 手动运行，Python 版可以批量改变参数、记录轨迹、导出解，并在同一脚本中计算多类指标。

### 5.2 支持增量加速与回归测试

项目后续在 `conop_py/incremental.py` 中开发了 `FastOrdinalState`，通过维护 per-section ordinal 缓存、事件位置数组和 numba JIT 内核，将 ordinal-only 模拟退火从十余秒级压缩到约一秒级。`tests/test_regression.py` 则锁定了 ordinal、level、eventual、随机扰动回滚和 schema 校验等行为，避免后续改动破坏已对齐的代价函数。

### 5.3 支持多目标扩展

`combined_misfit()` 将 Ordinal、Level 和 Eventual 通过权重组合为统一目标，例如 `3×ordinal + level + eventual`。这使项目可以从单一“最小化某个罚分”转向多目标权衡，讨论不同地层解释目标之间的 trade-off。若没有透明的 Python 代价函数，这类扩展很难在原版 CONOP9 上完成。

## 6 讨论

Python 复现工作最重要的成果，不是简单得到比 CONOP9 更低或更高的数值，而是建立了一套可解释的实验平台。在原版程序中，best fit 是最终输出；在复现版中，每一分 penalty 都能追溯到具体剖面、具体事件和具体算法规则。对于地层学而言，这种可追溯性有助于区分三类问题：真实的地层矛盾、观测数据不足导致的不确定性，以及算法或参数造成的搜索误差。

此外，复现结果也提醒我们，地层大数据分析不应只把计算程序当作黑箱。CONOP 的优化目标包含多层地质假设，例如 FAD/LAD 的先后、共存约束、锚点顺序、marker 固定方式和不同 misfit 的权重。每一个假设都会影响最终复合序列。因此，在课程论文和课堂汇报中，应同时展示算法结果和假设条件。

## 7 结论

1. 项目成功将 CONOP 输入文件解析为 Python 数据结构，并用 `(entity_id, event_type)` 表示复合序列中的 120 个事件。
2. `ConopContext` 统一保存事件位置、剖面观测和 taxon 延限信息，使 Ordinal、Level 与 Eventual 代价函数可共享预计算结构。
3. Python 版 Ordinal misfit 在原版解上得到 367，Level misfit 得到 237，均与 CONOP9 完全一致。
4. Eventual misfit 得到 335，较 CONOP9 的 353 低 18，主要可能源于 PAV 保序回归中偶数块中位数选择等实现细节。
5. 透明复现为后续模拟退火加速、多目标权衡、接受准则变体和解不确定性分析提供了可靠基础。

## 参考文献

Sadler P M, Cooper R A. 2003. Best-fit intervals and consensus sequences. In: Harries P J ed. High-Resolution Approaches in Stratigraphic Paleontology. Dordrecht: Springer, 49-94.

Sadler P M, Kemple W G, Kooser M A. 2003. CONOP9 programs for solving the stratigraphic correlation and seriation problems as constrained optimization. In: Harries P J ed. High-Resolution Approaches in Stratigraphic Paleontology. Dordrecht: Springer, 133-148.

Barlow R E, Bartholomew D J, Bremner J M, Brunk H D. 1972. Statistical Inference under Order Restrictions. New York: Wiley.

Best M J, Chakravarti N. 1990. Active set algorithms for isotonic regression. Mathematical Programming, 47: 425-439.

项目组. 2026. CONOP Python 复现代码[CP/OL]. geology-big-data-hw-main/conop_py.

项目组. 2026. Eventual misfit 溯源报告[Z]. geology-big-data-hw-main/results_py/eventual_diag/report.txt.

