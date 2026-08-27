# Paper 2 D104：last-passage 数学指导下的 CHL endpoint 分解（结果前冻结）

日期：2026-08-25  
状态：`PREREGISTERED_BEFORE_D104_CHL_OUTCOMES`

## 1. 数学问题

令径向距离为 `X_t`，running maximum 为 `M_t`，drawdown 为 `D_t=M_t-X_t`。冻结 RD 事件中：

- `tau_delta=inf{t:D_t>=delta}` 是可在线确认的首次 drawdown stopping time；
- `rho_delta=sup{s<=tau_delta:X_s=M_s}` 是 `tau_delta` 前最后一次径向最大值，只能由后续 drawdown 回标，不是 stopping time。

既有 CHL 结果是 `rho_delta` 相对已完成 `start→rho` 路径的事后秩富集。D102/D103 已经表明，简单 online CHL rank/background 不能跨 Goto 与 RFBO 稳定增加 renewal hazard。本轮不再重复 hazard 扫描，而回答两个尚未正式拆开的数学问题：

1. endpoint 双尾来自“所有径向 record 时刻本来就偏向 CHL 极端”，还是来自“未来 drawdown 从这些 records 中选择了最后一个 record”？
2. CHL 极端是定位在 last maximum `rho`，还是也延续到首次 drawdown trigger `tau`？

## 2. 数据、尺度和不重跑边界

- Goto：正式 Paper 2 主尺度 RD100，复用 canonical `start/endpoint/pretrigger/trigger` 索引、正式 strict-middle CHL 事件表和已冻结 4-km continuous-dense CHL cache；
- RFBO：复用 D27 movement-only RD 与 D28 CHL，尺度固定为 250/500/1000/2000 m；不加入 4000 m，不选择新尺度；
- Uesaka：只作结果盲分辨率资格审计。RD100 中虽然 3,674 个非 terminal 事件的 `endpoint→trigger` 中位为 10 s，但 daily 4-km observation cell 在 `rho→tau^-` 间仅 4.30% 改变，不能承担正式 CHL 时序检验；不打开其 CHL 效应；
- 不重跑 RD，不改变 delta，不下载新数据，不访问网络，不读写 Paper 1，不修改 Paper 2 正文。

Goto/RFBO 的 movement geometry 仅用于在既有冻结事件内部重建严格径向 record fixes；事件起止、endpoint 和 trigger 均以冻结表为准。若重建的最终 strict record 不等于冻结 endpoint，事件标为结构失败并退出，不修正事件。

## 3. CHL observation token 与 tail 定义

- 沿用正式 4-km daily cell、连续相同 cell 压缩为一个 run、run 内 CHL 聚合及 tie-aware average rank；
- high：rank fraction `>0.8`；low：`<=0.2`；union：high 或 low；
- 至少 3 个 finite CHL runs 且存在变化；
- Goto沿用 continuous-dense bilinear cache；RFBO沿用 D28 GPS-fix daily CHL，不更换采样核。

## 4. Family A：endpoint excess 的精确 last-passage 分解

对每个 `start→rho` step，以最终完整 step 的 tail flags 定义：

- `E = I(rho in tail) - mean_all_runs I(run in tail)`：既有 endpoint excess；
- `R = mean_record_runs I(run in tail) - mean_all_runs I(run in tail)`：径向 record-sampling 分量；
- `L = I(rho in tail) - mean_record_runs I(run in tail)`：last-record selection 分量。

逐事件必须数值满足 `E=R+L`。strict radial record 从 origin 后第一个 fix 开始，定义为相对 origin 的径向距离严格刷新 running maximum；同一 CHL run 内多个 record fixes 只计一个 record run。正式分解要求至少 2 个 distinct record runs；不足者只进入 E 复现审计，不进入 R/L 推断。

high、low、union 三项全部报告。Family A 的解释仅是事后 event-selection decomposition：

- `R>0` 表示径向进展本身较常经过该 step 的 CHL 极端；
- `L>0` 表示最终被未来 drawdown 回标的 last record 比同一 event 的其他 record runs 更偏极端；
- 两者都不能说明鸟直接感知 CHL。

## 5. Family B：last maximum `rho` 与 stopping trigger `tau` 的定位

对每个 `start→tau` 完整确认窗口，用同一窗口的 tail flags 计算：

- `E_rho_tauwindow = I(rho in tail)-mean_all_runs`；
- `E_tau_tauwindow = I(tau in tail)-mean_all_runs`；
- `C = I(rho in tail)-I(tau in tail)`，并验证 `C=E_rho_tauwindow-E_tau_tauwindow`。

主分析保留 rho/tau 落在同一 observation cell 的事件（它们提供真实的零对比）；另报告 distinct-cell 覆盖和同样的描述统计，但 distinct-cell 子集不替代主门。Family B 是 event-conditioned localization，不冒充 online hazard。严格可预测的 risk-intensity 结论继续引用 D102/D103，不另扫窗口、rank 门或模型。

## 6. 推断、null 与多重校正

- 观察效应先在 Goto track / RFBO BandNum 内平均，再 unit 等权；
- 95% CI：固定种子 5,000 次 unit cluster bootstrap；
- strong null：沿用 segment-common nonzero circular phase shift；movement、event boundary、record-run positions、rho/tau positions 和每个窗口长度固定，整段 CHL token 顺序共同错位后重新计算 tail 与各分量；199 次正式 phase replicates；
- Goto：Family A 的 2 分量×3 tails 共 6 项 Holm；Family B 的 3 estimands×3 tails 共 9 项 Holm；
- RFBO：相应 family 再乘 4 个冻结尺度后分别 Holm；
- `E`复现值与覆盖统计是结构审计，不参与追阳性判据。

单项支持同时要求：unit-bootstrap 95% CI 位于预期正方向，且 family 内 Holm phase `p<=.05`。负方向完整报告，不改成单侧反向假说。

## 7. 冻结裁决

Family A 对每个 tail 分别裁决：

- Goto与至少两个相邻RFBO尺度的 `R` 通过、`L`不形成同样复制：`RECORD_SAMPLING_DOMINANT`；
- `L`通过、`R`不形成同样复制：`LAST_RECORD_SELECTION_DOMINANT`；
- 两者都跨系统复制：`BOTH_COMPONENTS_REPLICATED`；
- 均不满足：`NO_REPLICATED_COMPONENT_LOCALIZATION`。

Family B 只有 `C=rho-tau` 在 Goto 与至少两个相邻 RFBO 尺度通过，才称 `LAST_MAXIMUM_LOCALIZATION_REPLICATED`；否则只能称数据集/尺度特异或未复制。

无论结果如何，D104都不得改写成：CHL 是感觉 cue、CHL 等于食物、rho 是鸟当时已知的 endpoint、tau 是捕食行为、已找到 Lévy 生成机制、或已推翻 Lévy 机制。

## 8. 停止规则

- smoke 只核对 schema、endpoint 复现、`E=R+L`、`C=E_rho-E_tau`、null 非退化、CPU affinity；
- 正式结果打开后不新增 tail cutoff、delta、物种、lag、天气、绝对 CHL 交互或行为分层；
- 若某系统无法复现冻结 endpoint / CHL token 规则，保留失败审计并从跨系统裁决中退出，不修改定义救结果。
