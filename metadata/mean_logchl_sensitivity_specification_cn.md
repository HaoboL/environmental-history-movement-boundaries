# D130 路径内 mean(log CHL) 最小敏感性：结果盲预注册

日期：2026-08-27

## 决策问题

检验 Figure 3 的同事件双参照结论是否依赖于以路径内部 `median(log CHL)` 表示 bout-level CHL background。主定义保持中位数；本分析只作一个预先限定的 `mean(log CHL)` 替代。

## 冻结输入与唯一改动

- 事件、物种、尺度和纳入规则与 D127 完全相同：Goto Crozet RD100，以及 Laysan 0.5、1、2 km。
- endpoint exclusion、连续原生 CHL cell-run 去重复、animal/track cluster 单位、`L_low`、`L_high`、`L_union`、step length 和事件连接均不改变。
- 只把 `interior_median_logchl` 替换为事件表中已经存在的 `interior_mean_logchl`。这里的 mean 是相同 endpoint-excluded interior runs 上的 `mean(log CHL)`；不加入空间 buffer、时间或距离权重，也不使用 raw-CHL arithmetic mean。
- 联合模型继续在 animal 内标准化 absolute background 与 log length，并同时纳入 animal-demeaned `L_low` 与 `L_high`。
- 不重跑 RD、CHL 下载、事件构建、Laysan last-record reconstruction 或 phase null。

## 推断与门槛

- 使用固定种子 1300827 和 20,000 次 animal-cluster bootstrap。
- 每个系统输出 absolute coefficient、`L_low`、`L_high`、标准化 `L_union` 及 95% CI，并输出与 D127 相同的 3×3 absolute-background × length 条件格。
- 总门沿用 PR-V3-017：
  1. Goto 与至少两个相邻 Laysan 尺度的 absolute coefficient 95% CI 上界小于 0；
  2. 同一组系统的 conditional `L_union` 95% CI 下界大于 0，且每个系统至少 6/9 条件格为正；
  3. 通过系统中的 `L_low` 点估计方向为正，并至少 Goto 与一个 Laysan 尺度的 `L_low` 95% CI 下界大于 0。

## 结果使用和停止规则

- 若通过：只作为主定义的最小敏感性，在 Results 或 SI 用一句话和一张紧凑补充表报告，不升级为独立发现。
- 若不通过：主结论收缩为 median-defined bout background，并如实报告替代结果。
- 完成本次分析后停止扩展 background summary；不得在读取结果后增加新门槛或继续尝试其他汇总量。

## 授权与状态

用户于 2026-08-27 明确要求执行此前全部待处理问题，视为对 PR-V3-017 冻结方案的确认。状态：`PREREGISTERED_BEFORE_RESULT_READ__READY_TO_RUN`。
