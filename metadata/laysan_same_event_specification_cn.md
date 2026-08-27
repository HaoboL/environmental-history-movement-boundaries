# D127 同事件双参照整合：结果盲预注册

日期：2026-08-26  
状态：`PREREGISTERED_BEFORE_JOINT_EFFECT_READOUT`

## 问题

在同一冻结 RD movement event 内，绝对 CHL 背景与回顾性的 last-passage 分量是否提供不可互换的信息：

1. `absolute-rich--short` 在加入 last-passage 分量后是否保留；
2. last-passage 的 `L_union = endpoint union-tail - strict radial-record union-tail mean` 是否并非只由绝对背景或 step 长度的样本构成造成。

这里的“整合”只指同事件统计非冗余，不解释为鸟在线计算 CHL 秩、不解释为感觉线索，也不解释为因果。

## 冻结输入与队列

- Goto：RD100；D104 `family_a_event_metrics` 与 D43 `absolute_chl_real_step_event_table` 按冻结 `event_id` 一一连接。
- Laysan：只使用 D42/D43 事先通过的相邻尺度 500、1000、2000 m；从既有 start-to-endpoint 路径点和 shift=0 CHL cell-run 表重建 strict radial records，不重跑 RD、不重选尺度。
- 主背景量：排除 endpoint 的 `interior_median_logchl`；主长度量：冻结 `step_length_km`。
- 推断单位：animal；Goto 使用 `cluster`，Laysan 使用 band/individual。

结果盲 key-only 审计已经完成：Goto 11,558 个连接事件；Laysan 500/1000/2000 m 分别 4,242/3,573/2,757 个 CHL 与路径均可连接事件。该审计未读取新的联合效应。

## Laysan last-passage 重建

对每个冻结 start-to-endpoint 路径：

1. 在球面单位向量的弦长坐标中，从 start 向后扫描严格增加的径向纪录；
2. 要求冻结 endpoint 是最后一个 strict radial record；
3. 将 strict record fixes 映射到既有 4-km daily CHL cell runs；同一 cell run 只计一次；
4. 至少 3 个 CHL runs、至少 2 个不同 strict-record CHL runs、CHL 非全相等才进入分解；
5. 按 D104 完全相同的 mid-rank 规则定义 high `>0.8`、low `<=0.2`、union，并计算 `E=R+L`。

## 分析 A：两层量在同一长度模型中是否非冗余

每个 dataset×scale 内，先在 animal 内中心化/标准化绝对背景；模型为：

`z_animal(log step length) ~ z_animal(absolute background) + L_low + L_high`

系数以 animal cluster bootstrap（20,000 次）给 95% percentile CI。另报告用 `L_union` 代替两个单尾的稳定性模型，但不用它替换预注册主模型。

主方向：

- absolute background 系数 `<0`；
- `L_low` 系数 `>0`（局地 relative-low 与较长 movement 对应）；
- 不要求 `L_high` 有固定方向。

## 分析 B：last-passage excess 是否跨绝对背景与长度条件存在

为避免连续协变量函数形式决定结论，在每个 animal 内分别把 absolute background 与 log length 冻结为 tertiles，形成最多 9 个同事件格。对每个 animal 先在每个非空格求 `L_union`，再对该 animal 的非空格等权，最后对 animal 等权。以 animal cluster bootstrap（20,000 次）给 CI，并报告 3×3 每格方向，防止总体值由单一 rich/poor 或 long/short 格驱动。

这一步的阳性门是标准化后的 `L_union` 95% CI 下界 `>0`；同时至少 6/9 个可估格方向为正。结构 phase null 的排伪由 D104（Goto）和本轮 Laysan 同定义 segment-common circular phase 检验承担；条件化结果本身不重新发明 null。

## 正式成功门

统一双参照现象只有同时满足才通过：

1. Goto：条件标准化 `L_union` CI 下界 `>0`，且 absolute background 长度系数 CI 上界 `<0`；
2. Laysan：至少两个相邻预冻结尺度同时满足上述两项；
3. Goto 与通过的 Laysan 尺度中，`L_low` 长度系数方向为正；跨系统不要求每一个 CI 都排除 0，但必须至少 Goto 和一个 Laysan 尺度 CI 下界 `>0`；
4. Laysan 新重建的 `L_union` 在至少两个相邻尺度同时通过 animal bootstrap 正 CI 与同尺度 family-wise Holm 校正后的 segment-common phase `p<=.05`；
5. `E=R+L` 最大绝对误差 `<=1e-12`，路径 endpoint 必须是最后 strict record，输入/事件覆盖/CPU affinity 全部审计通过。

失败时只能写两项统计在不同分析中并列，不能称同事件双参照组织。

## 禁止事项

- 不改 RD、CHL、tail、尺度或 animal 定义；
- 不新增天气、lag、感觉 proxy 或物种；
- 不修改 Paper 1；
- 不使用 CPU1；
- 不把结果称为 online trigger、感觉机制、捕食成功、Lévy 生成或因果。

