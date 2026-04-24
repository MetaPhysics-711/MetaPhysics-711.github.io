---
title: "误差分析"
abbrlink: 3276ef82
date: 2025-03-29 15:02:36
tags: "计算物理"
categories:
  - Physics Stack
  - Computational Physics
top_img: https://res.cloudinary.com/digumuwth/image/upload/top6_p7qqd5?_a=BAMAJaWO0
cover: https://res.cloudinary.com/digumuwth/image/upload/cover7_otwwlk?_a=BAMAJaWO0
---
# 问题引入
考虑一个多步运算的计算流程：
```mermaid
    Input --> U1 --> U2 --> ... --> Un --> Output
```
假设每个逻辑步骤 $U_i$ 相互独立，但每一步出错概率都是$p$，那么系统输出完全正确的联合概率为：
$$ P_{\text{correct}} = (1 - p)^n $$
由概率论数列收敛性知识可知，$p \in (0, 1)$，$(1 - p)<1$，$(1 - p)^n \rightarrow 0, \;  \text{as} \; n\rightarrow \infty$。

实际上，$n = 1000$，$p = 0.0001$，$P≈0.9048$。这已经不是很理想了。

# 误差的种类和来源
# 绝对和相对误差
# 误差的传播和估计
# 问题举例