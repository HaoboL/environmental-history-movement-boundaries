# D128 边界放置的宏观后果反事实：结果盲预注册

日期：2026-08-26  
状态：`PREREGISTERED_BEFORE_COUNTERFACTUAL_READOUT`

> **2026-08-26结构自检修订（效应读出前）：** 原稿把“每个冻结事件的 `tau` 按事件顺序连接”误当成合法 partition。只读索引审计发现，正式宏观尾部分布宇宙实际是 129 个 segment（233 是上游 complete-case risk-table 数，只有 129 个同时满足原 `MIN_BLOCKS/MIN_STEPS` 尾拟合门）；其中 70/129 个 segment 至少有一对相邻事件的 `tau` 非递增，共 124 对，52 个 segment 还存在重复 `tau`。原因是从前一 `rho` 开始的下一 RD event 可在前一 event 完成确认之前被确认。因此 `tau` 是每个回标事件的 confirmation mark，不是一套可直接替换 `rho` 的 renewal boundaries；强行连接会产生倒序或零长度。这是定义错误，不是结果阴性。正式主反事实据此限定为唯一保持事件数、时间顺序和正长度的 `rho` vs uniformly sampled eligible strict-record boundary`。`tau` 只报告 start→tau 的重叠 confirmation-reach 诊断，明确不得拟作 renewal/Brownian/Lévy partition，也不进入成功门。下文所有“233”改为“129”，所有 tau partition/attenuation 门由这一修订替代。

## 问题

Goto 中已经存在的 `rich--short / rich--Brownian-like, poor--Lévy-like` 经验关联，是否对 movement boundary 的定义敏感；具体地说，未来 drawdown 回标的 last maximum `rho` 是否比 online confirmation point `tau` 或普通 eligible radial record 更能产生该宏观表型。

这检验的是统计边界放置的后果，不声称动物执行 `rho`，也不改动 Paper 1 的 RD 定义或方法论文结论。

## 冻结宇宙

- 仅 Goto、RD100；不重跑 RD。
- 只使用既有 TWO_SIDED_COX 正式 `segment_rate_tail_metrics.csv` 中 129 个通过原 `MIN_BLOCKS/MIN_STEPS` 门的 segment；segment mean logCHL 和 segment/track 身份冻结。
- 每个非终止事件保留同一事件数和顺序，构造三套边界：
  1. `rho`：冻结 RD endpoint；
  2. `record`：每个事件从 start→rho 的全部 strict radial records 中等概率抽取一个；999 个 segment-preserving 反事实重复；
  3. `tau` 仅保留为每个事件 start→tau confirmation reach 的结构诊断，不构成替代分割。
- `rho` 与每次 `record` 边界按时间顺序形成完整分割：第一个长度为该 segment 冻结起点到首边界，后续为相邻边界间球面距离。不得把互相重叠的 start→tau 半径冒充新的 partition。

## 锚定与结构审计

- `rho` 重建长度必须逐事件复现冻结 `step_length_km`（最大绝对差预设 `<0.01 km`；超过则停止并诊断）。
- `rho` 重新拟合的 segment `AIC_exp-AIC_lomax` 与原正式表最大绝对差 `<1e-6`，且 CHL--support Spearman 复现到 `<1e-10`。
- 每次 `record` 边界必须严格按时间递增、在原 segment 内、事件数与 `rho` 相同；零/非法距离需审计，不能静默删除后改变 segment 宇宙。
- `tau` 必须审计事件顺序交叉数、重复数和 start→tau reach 的正值/覆盖；不得再要求它形成不存在的 renewal partition。

## 宏观统计

`rho` 与每次 eligible-record 分割在相同 129 segment 中用原正式函数拟合 exponential 与 Lomax，定义：

`Lomax support = AIC_exp - AIC_lomax`

并计算：

1. `rho(mean segment logCHL, Lomax support)`；预期为负；
2. `rho(mean segment logCHL, median log step length)`；预期为负（rich--short）。

对 `rho` 使用 track cluster bootstrap 20,000 次。对 999 个 `record` 反事实，报告统计量分布、2.5/50/97.5百分位和真实 `rho` 在该分布中的随机化 p 值。`tau` 的重叠 confirmation reach 只报告描述性 CHL 相关，不进入 Brownian/Lévy 推断。

因为 `rho` 与 `record` 每个 segment 的事件数完全相同，renewal density 在构造上不变；本实验不能、也不会把 `rich--dense` 的强弱归因于边界位置。该不变量必须在结果中明确报告，禁止暗示检验了不可能改变的量。

## 正式成功门

只有同时满足以下条件，才能称“经典宏观 phenotype 的 tail/length 部分由 environment-linked last-passage boundary placement 实质增强”：

1. `rho` 的 CHL--Lomax support 关联方向为负且 track-bootstrap CI 上界 `<0`；
2. `rho` 的 CHL--median-length 关联方向为负且 CI 上界 `<0`；
3. 真实 `rho` 的两项统计均比 eligible-record 反事实中位数更负，且至少一项落在反事实分布下 2.5% 尾（单侧随机化 `p<=.025`）；
4. 所有锚定、record时间/事件数审计和tau交叉/重复审计通过。

若只满足 `rho` 关联而 record 不削弱，结论为宏观关联对 eligible-record 边界放置稳健，不能把它解释为 last-passage selection 的后果。若 record 更强，则明确否决该升级假说。`tau` 因不构成 partition，无论其描述性方向如何都不用于支持或否决该宏观门。

## 禁止事项

- 不改变 233 segment、RD100、tail模型、拟合函数、最小事件数或环境量；
- 不从 999 次结果中挑 record 规则；
- 不把统计反事实写成动物在线策略或因果干预；
- 不修改 Paper 1，不使用 CPU1。
