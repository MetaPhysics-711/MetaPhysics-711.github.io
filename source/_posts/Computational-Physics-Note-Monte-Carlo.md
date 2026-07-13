---
title: 蒙特卡洛(Monte Carlo)方法
categories:
  - Physics Stack
  - Computational Physics
mathjax: true
abbrlink: eeb6438a
date: 2025-03-27 23:58:05
tags: "计算物理"
top_img: https://res.cloudinary.com/digumuwth/image/upload/top12_wa15gc?_a=BAMAJaWO0
cover: https://res.cloudinary.com/digumuwth/image/upload/cover15_ypj8ze?_a=BAMAJaWO0
---

# 简介

蒙特卡罗方法也称统计模拟方法，是指使用随机数（或者更常见的伪随机数）来解决很多计算问题的方法。它的工作原理就是两件事：不断抽样、逐渐逼近。

# 相关数学基础

## 条件概率

$P(A|B)$：在随机事件$B$发生的条件下，随机事件$A$发生的概率。

全概率公式：$P(B) = \sum_{i=1}^n P(B|A_i) \cdot P(A_i)$

- 证明思路：$P(B|A_i) \cdot P(A_i)$等于$P(BA_i)$，再利用概率的可加性进行加和。

贝叶斯(Bayes)公式：$P(A_i|B) = \frac{P(A_i B)}{P(B)} = \frac{P(B_i|A)P(A)}{\sum^n_{j=1} P(B|A_j) \cdot P(A_j)}$

## 随机变量的特征值

### 期望

- 离散分布：$E(Y(X)) = Y(X) \cdot p\,(Y(X))$

- 连续分布：$E(Y(X)) = \int_{-\infty}^{\infty} y(x) \cdot f(x) \, dx$，$f(x)$：概率密度函数

### 方差 

$D(X) = E[(X - E(X))^2] = E(X^2) - E^2(X)$

- 连续分布：$D(Y(X)) = \int_{-\infty}^{\infty} (y(x) - E(Y))^2 \cdot f(x) \, dx = E(Y^2) - (E(Y))^2$，$f(x)$定义同上。

## 常见分布

### 二项分布 (Binomial Distribution)

单次实验有两种结果，发生目标事件（概率为$p$）或者不发生（概率为$1-p$）。那么$N$次实验发生$k$次目标事件的概率、期望、方差为：
$$P(N, k) = \frac{N!}{k!(N - k)!} p^k (1 - p)^{N - k} \\
E(k)    = Np, \quad D(k) = Np(1 - p)$$

- 对应物理实例：考虑一个由 $N$ 个全同自旋 $1/2$ 粒子组成的系统，每个粒子的磁矩为 $m$，外磁场强度为 $B$，自旋可以和磁场同向或反向，对应能量为$E = \mp mB$。若系统中有 $n$ 个粒子自旋同向，$N-n$ 个粒子自旋反向，则总能量为 $U(n) = -nE + (N-n)E = -(2n - N)E$。熵为 $\sigma(n) = \log \Omega(n) = \log \binom{N}{n} = \log \left( \frac{N!}{n!(N-n)!} \right)$。

### 泊松分布 (Poisson Distribution)

泊松分布是二项分布的一种特殊极限，当单次实验 $p \rightarrow$ 0 、$n \rightarrow \infty$，且 $Np\,$适中（有限）时，在相同时间内，随机过程发生$k$次的概率为：
$$
P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}, \quad k = 0, 1, 2, \ldots
$$
其中$\lambda > 0$ 为平均发生次数（速率参数）。

泊松分布的方差和期望均是$\lambda$。

- 对应的物理实例：放射性同位素 $^{137}\text{Cs}$ 的半衰期大约是 27 年，那么每秒钟衰变的可能性大约为 $D = \frac{1\,\text{s} \cdot \ln 2}{27\,\text{years}}\approx 8.14 \times 10^{-10} \, \text{s}^{-1}$。如果我们有 $1\,\mu\text{g}$ 的 $^{137}\text{Cs}$ 样品，大约有 $N = 10^{15}$ 个核子。那么此时每秒约有 $N D = 8.14 \times 10^{5}$ 次衰变。
此时我们可以认为该样品的衰变数统计符合泊松分布。

### 正态分布 (Normal Distribution)

自然科学中（同时也是生活中）最常见的分布形式？

高斯分布的概率密度函数（PDF）：
$$
f(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(x - \mu)^2}{2\sigma^2}}, \quad x \in \mathbb{R}
$$
正态分布的期望和方差分别为：
$$
E(X) = \mu, \quad D(X) = \sigma^2
$$
- 对应的物理实例：量子力学中的一维谐振子，哈密顿量为 $H = \frac{p^2}{2m} + \frac{1}{2} m\omega^2 x^2$。它的基态波函数 $\psi_0(x)$ 可以根据升降算符求得，根据 $a\psi_0(x) = 0$，
我们可以求解：
  $$
  \left( x + \frac{\hbar}{m\omega} \frac{d}{dx} \right) \psi_0(x) = 0
  $$
  由此，我们可以得到它符合均值为 0 的高斯基态波函数 $\psi_0(x) = Ce^{-\frac{m\omega x^2}{2\hbar}}$。

## 中心极限定理
中心极限定理指出，**无论原始随机变量的分布如何**，其样本均值的标准化形式在样本量 $n$ 足够大时，将近似服从标准正态分布。这是统计学中连接概率分布与正态分布的桥梁。

设 $X_1, X_2, \dots, X_n$ 是独立同分布（i.i.d.）的随机变量，且$E(X_i) = \mu$，$D(X_i) = \sigma^2 < \infty$。定义样本均值为：$\overline{X_n} = \frac{1}{n} \sum_{i=1}^n X_i$，则当 $n \to \infty$ 时，标准化随机变量：$Z_n = \frac{\overline{X_n} - \mu}{\sigma / \sqrt{n}}$的分布收敛于标准正态分布 $N(0,1)$，即：
$$
Z_n \xrightarrow{d} N(0,1)
$$

# 随机数产生器
（不确定是否重要……有空再来填这个坑吧）

# 随机抽样

## 离散随机变量
离散随机变量的抽样相对简单，我们可以根据各变量的频率及其概率，从 $[0,1)$ 均匀分布中分区间直接抽样。

具体方法：设随机变量为$X_1, X_2, ..., X_n$，相应概率为$P_1, P_2, ..., P_n$。

定义累积概率 $\xi$：$\xi_0 = 0$，$\xi_i = \sum_{j=1}^i p_j, i = 1, 2, ..., n$。接着我们从$[0, 1)$中均匀分布抽样得到$x^*$，如果满足：
$$\xi_{k-1} \leq x^* \leq \xi_k$$
那么此时抽样的随机变量为$X_k$。

辅助理解的图例（谢谢伟大的陈海老师的ppt！）：![](https://res.cloudinary.com/digumuwth/image/upload/%E5%B1%8F%E5%B9%95%E6%88%AA%E5%9B%BE_2025-07-15_182731_aopvfz?_a=BAMAJaWO0)

## 连续随机变量

### 直接抽样法
### 变换抽样法
### 舍选抽样法

# 蒙特卡洛方法的具体应用

## 积分运算

### 一维积分运算

#### 投点法
#### 平均值法
#### 重要抽样法

### 多维积分运算

## 统计力学中的运用

### Metropolis-Hastings算法

## 量子力学中的运用

### 路径积分

## 模拟退火算法

## 随机行走的生长模拟