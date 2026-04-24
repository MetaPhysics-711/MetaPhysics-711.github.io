---
title: 经典数值计算
abbrlink: 994c071a
date: 2025-04-07 09:07:09
tags: "计算物理"
categories: 
  - Physics Stack
  - Computational Physics
top_img: https://res.cloudinary.com/digumuwth/image/upload/top7_uxngzh?_a=BAMAJaWO0
cover: https://res.cloudinary.com/digumuwth/image/upload/cover1_w3xlh0?_a=BAMAJaWO0
---

# 插值问题

物理问题 $\rightarrow$ 离散化 $\rightarrow$ 计算机处理

离散的物理数据，但不知道实际函数$\rightarrow$ 构造近似函数作为实际函数$f(x)$的逼近

---

**插值问题的数学表述：** $f(x)$是定义在区间$[a, b]$的函数，$(x_0, y_0),\ (x_1, y_1),\ (x_2, y_2)...,\ (x_n, y_n)$是该区间上$n+1$个不同的点。我们需要构造一个便于计算的函数$P(x)$，s.t. 
$$
P(x_i) = y_i \tag{1.1}
$$
则称$P(x)$是$f(x)$的插值函数，$[a, b]$是插值区间，$(x_i, y_i)$是插值节点。同时在其他$x \neq x_i$点上，$P(x)$可以作为$f(x)$的近似。

## 多项式插值
设$P_n(x) = a_0 + a_1x + a_2x^2 + ... + a_nx^n,\;\;  a_0, a_1, ..., a_n \in \mathbb{R}$，且满足插值条件${1.1}$。

将插值条件写成线性方程组的形式：
$$
\begin{cases}
a_0 + a_1 x_0 + a_2 x_0^2 + \cdots + a_n x_0^n = y_0 \\
a_0 + a_1 x_1 + a_2 x_1^2 + \cdots + a_n x_1^n = y_1 \\
\vdots \\
a_0 + a_1 x_n + a_2 x_n^2 + \cdots + a_n x_n^n = y_n
\end{cases}
$$
这样可用矩阵语言表述：
$$
\begin{pmatrix}
y_0 \\
y_1 \\
\vdots \\
y_n
\end{pmatrix}
=
\begin{pmatrix}
1 & x_0 & \cdots & x_0^n \\
1 & x_1 & \cdots & x_1^n \\
1 & x_2 & \cdots & x_2^n \\
1 & x_n & \cdots & x_n^n
\end{pmatrix}
\begin{pmatrix}
a_0 \\
a_1 \\
\vdots \\
a_n
\end{pmatrix} 
$$

记作$(1.2)$。上述系数矩阵是方阵（这很显然），记作$X$。由线性代数的知识可知，当系数行列式$|X| \neq 0$时，齐次方程组$Xt = 0$只有零解，对应的非齐次方程组$Xt = b$只存在一组特解。

求解：
$$
Y = Xa\\
\Rightarrow a = X^{-1}Y
$$
其中$Y$和$a$分别表示$y_i$和$a_i$对应的的列向量，$X^{-1}$表示$X$的逆矩阵。


注意到系数矩阵恰巧时范德蒙德（Vandermonde）矩阵。由线性代数知识可知，Vandermonde行列式的值为（将对应矩阵记作V）：
$$
\det(V) = \begin{vmatrix}
1 & x_0 & \cdots & x_0^n \\
1 & x_1 & \cdots & x_1^n \\
\vdots & \vdots & \ddots & \vdots \\
1 & x_n & \cdots & x_n^n
\end{vmatrix}
= \prod_{0 \leq j < i \leq n}(x_i - x_j) \neq 0
$$

{% note info simple %}

Vandermonde行列式求解步骤：

已知$\det(V) = \begin{vmatrix}
1 & x_0 & \cdots & x_0^n \\
1 & x_1 & \cdots & x_1^n \\
\vdots & \vdots & \ddots & \vdots \\
1 & x_n & \cdots & x_n^n
\end{vmatrix}$

由行列式的性质：转置以及某行的$k$倍加到另一行，行列式的值不变。

因此做变换：先转置，再从最后一行开始$r_i - x_0r_{i-1}$（$r_i$表示$x_i$所在的行）得到

$$
\begin{align*}
\det(V^T) &= \begin{vmatrix}
1 & 1 & \cdots & 1 \\
x_0 & x_1 & \cdots & x_n\\
x_0^2 & x_1^2 & \cdots & x_n^2\\
\vdots & \vdots & \ddots & \vdots\\
x_0^n & x_1^n & \cdots & x_n^n
\end{vmatrix}\\
&= \begin{vmatrix}
1 & 1 & \cdots & 1 \\
0 & x_1-x_0 & \cdots & x_n-x_0\\
0 & x_1^2-x_1x_0 & \cdots & x_n^2-x_nx_0\\
\vdots & \vdots & \ddots & \vdots\\
0 & x_1^n-x_1^{n-1}x_0 & \cdots & x_n^n-x_n^{n-1}x_0
\end{vmatrix}\\
&= \begin{vmatrix}
x_1-x_0 & x_2-x_0 & \cdots & x_n-x_0\\
x_1(x_1-x_0) & x_2(x_2-x_0) & \cdots & x_n(x_n-x_0)\\
\vdots & \vdots &\ddots & \vdots\\
x_1^{n-1}(x_1-x_0) & x_2^{n-1}(x_2-x_0) & \cdots & x_n^{n-1}(x_n-x_0)
\end{vmatrix}\\
& = \prod_{i=1}^n(x_i-x_0)
\begin{vmatrix}
1 & 1 & \cdots & 1\\
x_1 & x_2 & \cdots & x_n\\
x_1^2 & x_2^2 & \cdots & x_n^2\\
\vdots & \vdots & \ddots & \vdots\\
x_1^{n-1} & x_2^{n-1} & \cdots & x_n^{n-1}
\end{vmatrix}
\end{align*}
$$

可以看到又出现了Vandermonde行列式，不过此时阶数$-1$。再进行一轮上述操作，又得到因子$\prod_{j=2}^n (x_j-x_1)$。以此类推，最终可得：
$$
\begin{align*}
\det(V^T) &= \det(V) = \prod_{i=1}^n(x_i-x_0)\prod_{j=2}^n (x_j-x_1)\prod \cdots (x_n-x_{n-1})\\
&=\prod_{0 \leq j < i \leq n}(x_i - x_j)
\end{align*}
$$

{% endnote %}

{% note flat %}

定理：对于$n+1$个互不相同的插值节点，满足插值条件的$n$次多项式插值函数存在且唯一。

{% endnote %}

### 多项式插值方法一 —— 拉格朗日插值法
#### 引入：线性插值
给定两个点$(x_0, y_0)$，$(x_1, y_1)$求解
$$
\begin{cases}
a_0+ a_1x_0 = y_0\\
a_0+ a_1x_1 = y_1
\end{cases}
$$
解出：
$$
\begin{cases}
a_0 = \frac{y_0x_1-y_1x_0}{x_1-x_0}\\
a_1 = \frac{y_0 - y_1}{x_0 - x_1}
\end{cases}
$$
可得
$$
P_1(x) = a_0+ a_1x = \frac{y_1x_0-y_0x_1}{x_0-x_1} + \frac{y_0 - y_1}{x_0 - x_1}x =\frac{x-x_1}{x_0 - x_1}y_0 - \frac{x-x_0}{x_0-x_1}y_1
$$
记 $l_0 = \frac{x-x_1}{x_0-x_1}, \; l_1 = \frac{x-x_0}{x_1-x_0}$，那么 $P_1 = l_0(x) y_0 + l_1(x) y_1 \; \Rightarrow \; l_0(x), l_1(x)$为插值基函数。
#### 推广：拉格朗日插值

{% note default simple %}

首先需要明白插值基函数的性质。它们满足在插值节点$x_k$上等于1，其余地方为0。即：
$$
l_k(x_i) = 
\begin{cases}
1 \quad i = k\\
0 \quad i\neq k
\end{cases}
$$

{% endnote %}

构造：
$$
  P_n(x) = \sum_{k=0}^n l_k(x) \cdot y_k, \quad l_k(x) = \prod_{j \neq k} \frac{x - x_j}{x_k - x_j}
$$
其中$l_k(x)$是基函数。

易验证拉格朗日插值基函数满足以下性质：
- 归一性：$l_k(x_i) = 
\begin{cases}
1 \quad i = k\\
0 \quad o/w
\end{cases}$
- 多项式次数：每个$l_k$是$n$次多项式。

#### 举例理解
上面的线性插值已经是一个例子。现在看$n = 2$：

对于$(x_0,x_1,x_2)$，
$$
\begin{aligned}
l_0(x) &= \frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)} \\
l_1(x) &= \frac{(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)} \\
l_2(x) &= \frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}
\end{aligned}
$$
更高阶同理。

#### 线性插值误差
{% note flat %}

定理：设 $f(x)$ 一阶连续可导，且 $f''(x)$ 在区间 $(a, b)$ 上存在，则对任意 $x \in [a, b]$，存在 $\xi \in [a, b]$ 使得：
$$
R(x) = f(x) - P_1(x) = \frac{f''(\xi)}{2}(x-x_0)(x-x_1)
$$
这里的 $P_1(x)$ 是线性插值的基函数。

{% endnote %}

证明：

$\because f(x_0) = P_1(x_0); \quad f(x_1) = P_1(x_1)$

$\therefore R(x) = f(x)-P_1(x)$存在两个零点$x_0, x_1$

$\therefore$ 可以设 $R(x) = k(x)(x-x_0)(x-x_1)$，且：
$$
f(x)-P_1(x) = k(x)(x-x_0)(x-x_1) \tag{1.2}
$$

构造辅助函数：
   $$\phi(t) = f(t) - P_1(t) - k(x)(t-x_0)(t-x_1)$$
则可以将${1.2}$看成$\phi(t)$在$t=x$时的取值，且$\phi(x) = 0$

则由以上讨论，$\phi(t)$ 在 $x_0,x,x_1$ 处有零点。由Rolle定理，$\phi'(t)$在 $(x_0, x_1),\, (x_1, x_2)$ 上存在零点，分别记作 $(\xi_1,\xi_2)$。再由罗尔定理，$\phi''(t)$ 在 $(\xi_1,\xi_2)$ 有零点，记作 $\xi$。

$\therefore \phi''(\xi) = 0$ 对 $\phi(t)$ 两边求二阶导得：
$$
\phi''(t) = f''(t) - P_1''(t) - 2k(x) = 0
$$
又 $\because P_1(t)$是线性插值基函数，仅有一次项 $\therefore P_1''(x) = 0$ $\quad \Rightarrow$ $f''(\xi) = 2k(x)$

再代回原式即可得证。

**一般形式的误差估计：** 对于 $n$ 次插值多项式 $P_n(x)$：

$$
R_n(x) = f(x) - P_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)
$$

其中：
$$
\omega_{n+1}(x) = \prod_{i=0}^n (x-x_i)
$$
可以看出：节点越密集 ($(x-x_i)$ 越小)，误差越小。

很多情况下，$f(x)$ 未知，$f^{(n+1)}$也难以得知，因此可采用以下方法估计误差：

使用 $n+2$ 个数据点进行两次插值运算：前 $n+1$ 点构造 $P_n(x)$；后 $n+1$ 点构造 $\hat{P}_n(x)$。

两次的误差分别是：
$$
R_n(x) = f(x) - P_n(x) = \frac{f^{(n+1)}(\xi_1)}{(n+1)!}(x-x_0)(x-x_1)\cdots(x-x_n)\\
\hat{R}_n(x) = f(x) - \hat{P}_n(x) = \frac{f^{(n+1)}(\xi_2)}{(n+1)!}(x-x_1)(x-x_2)\cdots(x-x_{n+1})
$$
整理可得误差关系式：
   $$
   f(x) - P_n(x) \approx \frac{x-x_0}{x_0-x_{n+1}}(P_n(x)-\hat{P}_n(x))
   $$
{% note info simple %}

这里说说得到这个误差关系式的具体步骤：

已知
$$
R_n(x) = f(x) - P_n(x) = \frac{f^{(n+1)}(\xi_1)}{(n+1)!}(x-x_0)(x-x_1)\cdots(x-x_n)\\
\hat{R}_n(x) = f(x) - \hat{P}_n(x) = \frac{f^{(n+1)}(\xi_2)}{(n+1)!}(x-x_1)(x-x_2)\cdots(x-x_{n+1})
$$
我们**假设 $f^{(n+1)}(\xi_1) \approx f^{(n+1)}(\xi_2)$，即高阶导数在 $\xi_1$ 和 $\xi_2$ 处的值相近**。于是可以将两个误差表达式中的高阶导数项消去：
$$
\frac{R_n(x)}{\hat{R}_n(x)} \approx \frac{(x - x_0)(x - x_1) \cdots (x - x_n)}{(x - x_1)(x - x_2) \cdots (x - x_{n+1})} = \frac{x - x_0}{x - x_{n+1}}
$$
$$
\Rightarrow R_n(x) \approx \frac{x - x_0}{x - x_{n+1}} \hat{R}_n(x)
$$
由于$\hat{R}_n(x) = f(x) - \hat{P}_n(x)$，则代入可得：
$$
f(x) - P_n(x) \approx \frac{x - x_0}{x - x_{n+1}} (f(x) - \hat{P}_n(x))
$$
将 $f(x)$ 的项移到左边：
$$
f(x) \left(1 - \frac{x - x_0}{x - x_{n+1}}\right) \approx P_n(x) - \frac{x - x_0}{x - x_{n+1}} \hat{P}_n(x)
$$
$$
\Rightarrow f(x) \frac{x_0-x_{n+1}}{x-x_{n+1}} \approx P_n(x)-\frac{x - x_0}{x - x_{n+1}} \hat{P}_n(x)
$$
$$
\begin{aligned}
\Rightarrow f(x) &\approx \frac{x-x_{n+1}}{x_0-x_{n+1}}[P_n(x)-\frac{x - x_0}{x - x_{n+1}} \hat{P}_n(x)]\\
&=\frac{x-x_{n+1}}{x_0-x_{n+1}}P_n(x) - \frac{x-x_0}{x_0-x_{n+1}}\hat{P}_n(x)\\
&=\frac{x_0-x_{n+1}+x-x_0}{x_0-x_{n+1}}P_n(x) - \frac{x-x_0}{x_0-x_{n+1}} \hat{P}_n(x)\\
&=P_n(x) + \frac{x-x_0}{x_0-x_{n+1}}[P_n(x) - \hat{P}_n(x)]
\end{aligned}
$$
移项整理即可得到误差关系式。

{% endnote %}

#### 拉格朗日插值法的优点/缺点：
- 优点：
  - 显示公式：无需解线性方程组，直接构造多项式。
  - 理论直观
- 缺点：
  - 计算效率低：新增节点需重新计算所有基函数。

#### 代码实现
```python
def Pn_x(x_variables, x_data, y_data):
    """
    计算拉格朗日插值多项式在 x_variables 处的值

    参数:
    x_variables: float 或 array-like, 插值点的自变量值
    x_data: list 或 array-like, 已知数据点的自变量值
    y_data: list 或 array-like, 已知数据点的因变量值

    返回:
    float 或 array-like, 插值多项式在 x_variables 处的值
    """
    n = len(x_data)  # 数据点数量
    result = 0.0  # 初始化结果

    for i in range(n):
        term = y_data[i]  # 初始化当前项
        # 计算拉格朗日基函数
        for j in range(n):
            if j != i:
                term *= (x_variables - x_data[j]) / (x_data[i] - x_data[j])
        result += term
    return result
```

### 多项式插值方法二——牛顿插值法
#### 引入
目的：避开拉格朗日插值法“新增节点需重新计算所有基函数”这一问题。

由数学知识可知任何一个多项式都可以表示为：
$$
N_n(x) = a_0 + \sum_{k=1}^n a_k \prod_{i=0}^{k-1} (x - x_i)
$$
在牛顿形式下，我们可以认为此时 $a_k$ 对应的基函数是 $k$ 的多项式$\prod_{i=0}^{k-1} (x - x_i)$，它们相互独立。因此有更高阶数加入的时候只需计算相应新的多项式系数即可。
#### 计算步骤
1. 将插值条件$P(x_i) = y_i$写成线性方程组的形式：
$$
\begin{cases}
y_0 = a_0 = f(x_0)\\
y_1 = y_0 +a_1 (x_1-x_0) = f(x_1)\\
y_2 = y_0 + a_1 (x_2-x_0) + a_2(x_2-x_0)(x_2-x_1) = f(x_2)\\
...\\
y_n = y_0 + a_1 (x_n-x_0) + a_2(x_n-x_0)(x_n-x_1) + ... + a_n(x_n-x_0)(x_n-x_1)...(x_n-x_{n-1}) = f(x_n)
\end{cases}
$$

2. 引入差商的概念：
   由$(y_0, y_1)$我们引入一阶差商
$$
F(x_0, x_1) = \frac{f(x_1) - f(x_0)}{x_1 - x_0} = a_1
$$
   类似可得更高阶的差商公式：
$$
F(x_0, x_1,x_2) = \frac{F(x_2, x_1)-F(x_0, x_1)}{x_2-x_0}\\
\cdots \\
F(x_0,...,x_n) = \frac{F(x_1,...,x_n)-F(x_0,...,x_{n-1})}{x_n-x_0}
$$

根据以上公式我们可以改写差商形式的牛顿插值多项式：
$$
\begin{aligned}
P_n(x) &= a_0 + \sum_{k=1}^n a_k \cdot \prod_{i=0}^{k-1} (x - x_i) \\
&= F(x_0) + \sum_{k=1}^n F(x_0,...,x_k) \cdot \prod_{i=0}^{k-1}(x-x_i)
\end{aligned}
$$
对应关系：$a_k = F(x_0, ..., x_k)$

{% note info simple %}

这里以 $a_2$ 为例，详细展开系数的计算：

易知
$$
\begin{aligned}
a_2 &= \frac{f(x_2)-f(x_0)-F(x_0, x_1)(x_2-x_0)}{(x_2-x_0)(x_2-x_1)}\\
&= \frac{f(x_2)-f(x_0)-\frac{f(x_1) - f(x_0)}{x_1 - x_0}(x_2-x_0)}{(x_2-x_0)(x_2-x_1)}\\
&= \frac{(f(x_2)-f(x_0))(x_1-x_0)-(f(x_1)-f(x_0))(x_2-x_0)}{(x_2-x_0)(x_2-x_1)(x_1-x_0)} \\
&= \frac{f(x_2)-f(x_0)}{(x_2-x_0)(x_2-x_1)} - \frac{f(x_1) - f(x_0)}{(x_1-x_0)(x_2-x_1)}\\
\end{aligned} \tag{1.3}
$$

而更常见的形式为：
$$
a_2 = \frac{\frac{f(x_2) - f(x_1)}{x_2 - x_1}-\frac{f(x_1) - f(x_0)}{x_1 - x_0}}{x_2-x_0} \tag{1.4}
$$
下证$(1.3)$和$(1.4)$等价性：
$$
\begin{aligned}
(1.3) &= \frac{(f(x_2)-f(x_0))(x_1-x_0)-(f(x_1)-f(x_0))(x_2-x_0)}{(x_2-x_0)(x_2-x_1)(x_1-x_0)}\\
&= \frac{(x_1f(x_2)-x_0f(x_2)-x_1f(x_0)+x_0f(x_0))-(x_2f(x_1)-x_0f(x_1)-x_2f(x_0)+x_0f(x_0))}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}\\
&= \frac{x_1f(x_2)-x_0f(x_2)-x_1f(x_0)-x_2f(x_1)+x_0f(x_1)+x_2f(x_0)}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}\\
&= \frac{f(x_2)(x_1-x_0)+f(x_1)(x_0-x_2)+f(x_0)(x_2-x_1)}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}
\end{aligned}
$$
$$
\begin{aligned}
(1.4) &= \frac{f(x_2)-f(x_1)}{(x_2-x_1)(x_2-x_0)} - \frac{f(x_1)-f(x_0)}{(x_1-x_0)(x_2-x_0)}\\
&= \frac{(x_1-x_0)(f(x_2)-f(x_1))-(x_2-x_1)(f(x_1)-f(x_0))}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}\\
&= \frac{x_1f(x_2)-x_1f(x_1)-x_0f(x_2)+x_0f(x_1)-(x_2f(x_1)-x_2f(x_0)-x_1f(x_1)+x_1f(x_0))}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}\\
&= \frac{x_1f(x_2)-x_0f(x_2)+x_0f(x_1)-x_2f(x_1)+x_2f(x_0)-x_1f(x_0)}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}\\
&= \frac{f(x_2)(x_1-x_0) + f(x_1)(x_0-x_2) + f(x_0)(x_2-x_1)}{(x_2-x_1)(x_2-x_0)(x_1-x_0)}
\end{aligned}
$$
可见 $(1.3) \iff (1.4)$ 。

而由定义可得
$$
F(x_2, x_1) = \frac{f(x_2) - f(x_1)}{x_2 - x_1}
$$
所以根据$(1.4)$，我们可以把 $a_2$ 写成
$$
a_2 \triangleq F(x_0, x_1, x_2) = \frac{F(x_1, x_2) - F(x_0, x_1)}{x_2-x_0}
$$
这就是上面给出的形式。
（所以怎么直接由$(1.3)$推出$(1.4)$……求大神评论区指点。我自己也会想想办法的。）

{% endnote %}

#### 牛顿差商的性质
1. 线性叠加：
   利用数学归纳法可以证明：
   $$
   F(x_0, x_1, ..., x_m) = \sum^m_{j=0}\frac{f(x_j)}{(x_j-x_0)\cdots(x_j-x_{j-1})(x_j-x_{j+1})\cdots(x_j-x_m)}
   $$
2. 线性性：
若$h(x) = \alpha f(x) + \beta g(x)$，则：
$$
H(x_0,...,x_n) = \alpha F(x_0,...,x_n) + \beta G(x_0,...,x_n)
$$
3. 对称性
差商值与节点排列顺序无关：
$$
F(x_0,x_1,...,x_n) = F(x_{\sigma(0)},x_{\sigma(1)},...,x_{\sigma(n)})
$$
  其中$\sigma$是任意排列
  
4. 差商和导数关系：
若 $f(x)\in C^n[a,b]$，则存在 $\xi\in(a,b)$ 使得：
$$
F(x_0,...,x_n) = \frac{f^{(n)}(\xi)}{n!}
$$
5. 差商表可加性：
新增节点$x_{n+1}$时：
$$
F(x_0,...,x_{n+1}) = \frac{F(x_1,...,x_{n+1}) - F(x_0,...,x_n)}{x_{n+1}-x_0}
$$
其中$F(x_1,...,x_{n+1})$可由性质1算出。
#### 代码实现
```python
def newton_interpolation(x_points, y_points, x):
    """
    牛顿插值法实现
    
    参数:
        x_points: 已知点的x坐标列表
        y_points: 已知点的y坐标列表
        x: 要计算插值的x值
        
    返回:
        插值结果y值
    """
    # 检查输入数据是否有效
    if len(x_points) != len(y_points):
        raise ValueError("x_points和y_points的长度必须相同")
    if len(x_points) == 0:
        raise ValueError("输入数据不能为空")
    
    n = len(x_points)
    
    # 构造差商表
    divided_diff = [[0] * n for _ in range(n)]
    for i in range(n):
        divided_diff[i][0] = y_points[i]
    
    for j in range(1, n):
        for i in range(n - j):
            divided_diff[i][j] = (divided_diff[i+1][j-1] - divided_diff[i][j-1]) / (x_points[i+j] - x_points[i])
    
    # 计算插值结果
    result = divided_diff[0][0]
    product_term = 1
    
    for i in range(1, n):
        product_term *= (x - x_points[i-1])
        result += divided_diff[0][i] * product_term
    
    return result
```
### 龙格(Runge)现象
指在高次多项式插值中，随着插值节点数量的增加，插值结果在区间端点附近出现剧烈振荡的现象。这种现象由德国数学家Carl Runge在研究多项式插值时首次发现。
- 这说明使用高次多项式插值并不总能提高准确性。

**解决方法：多次线性插值** 
将插值区间划分为若干子区间，在每个子区间上用**线性多项式**进行插值的方法。给定节点$(x_0,y_0),(x_1,y_1),...,(x_n,y_n)$，我们构造的插值函数应该满足：
- $S(x_i) = y_i$ （通过所有节点）
- $S(x)$在每个$[x_i,x_{i+1}]$上是线性函数
- 整体函数连续（$C^0$连续）

具体实现步骤：

1. 找到待插值点$x$及其所在的区间$[x_k,x_{k+1}]$
2. 计算：$P_1(x) = \frac{x-x_{k+1}}{x_k - x_{k+1}}y_k + \frac{x-x_k}{x_{k+1} - x_k}y_{k+1}$（参照拉格朗日插值的形式）

缺陷：$C^1$不连续，无法保证插值节点的光滑性。

## 样条插值
样条插值（Spline Interpolation）是一种**分段低次多项式插值**方法，通过强制拼接点处的光滑条件来保证整体曲线的平滑性。（相当于多次线性插值的一种解决方案）

### 原理
插值函数需要满足：
- $S_i(x_i)=y_i$, $S_i(x_{i+1})=y_{i+1}$（插值条件）
- 连续性：
   - $C^0$: $S_{i-1}(x_i) = S_i(x_i)$
   - $C^1$: $S'_{i-1}(x_i) = S'_i(x_i)$ 
   - $C^2$: $S''_{i-1}(x_i) = S''_i(x_i)$

最常用：三次样条插值函数，即在每段$[x_i,x_{i+1}]$上构造：
$$
S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3
$$
对$n+1$个插值点，总共有$n$个区间，那么总参数有$4n$个。

约束条件（加起来一共$4n$个）：
- 插值条件：提供 $n+1$ 个
- 连续性：提供 $3(n-1)$ 个
- 边界条件：提供2个

边界条件还能分为以下三种情况：
- 第一类边界条件：一阶导数值指定，即$S'(x_0)=y'_0$，$S'(x_n)=y'_n$
- 第二类边界条件：二阶导数值指定，即$S''(x_0)=y''_0$，$S''(x_n)=y''_n$
  - 自然边界条件：$S''(x_0)=S''(x_n)=0$
- 周期性边界条件：$S'(x_0)=S'(x_n)$，$S''(x_0)=S''(x_n)$

### 求解步骤
以上面提到的三次样条插值为例，对于$x \in [x_k, x_{k+1})$：
$$
\begin{aligned}
S_k(x) &= a_k + b_k h + c_k h^2 + d_k h^3 \quad (h = x-x_k) \\
S'_k(x) &= b_k + 2c_k h + 3d_k h^2 \\
S''_k(x) &= 2c_k + 6d_k h
\end{aligned}
$$
应满足条件：
- $S_k(x_k) = y_k \Rightarrow a_k = y_k$（插值条件）
- $S_k(x_{k+1}) = y_{k+1} \Rightarrow a_k + h_k b_k + h_k^2 c_k + h_k^3 d_k = y_{k+1}$ （其中 $h_k = x_{k+1}-x_k$）（ $C^0$ 连续性条件）
- $b_k + 2h_k c_k + 3h_k^2 d_k = b_{k+1}$（ $C^1$ 连续性条件）
- $2c_k + 6h_k d_k = 2c_{k+1}$（ $C^2$ 连续性条件）

求解该方程组：
令 $m_k = 2c_k$，则
1. 由二阶导数连续条件：
   $$ d_k = \frac{m_{k+1}-m_k}{6h_k} $$
2. 代入$C^0$连续性条件得：
   $$ b_k = \frac{y_{k+1}-y_k}{h_k} - \frac{h_k}{2}m_k - \frac{h_k}{6}(m_{k+1}-m_k) $$
3. 将上式代入$C^1$连续条件，整理得到：
$$
h_k m_k + 2(h_k+h_{k+1})m_{k+1} + h_{k+1}m_{k+2} = 6\left(\frac{y_{k+2}-y_{k+1}}{h_{k+1}} - \frac{y_{k+1}-y_k}{h_k}\right)
$$
4. 加上两个边界条件的约束 $\Rightarrow$ 建立线性方程组求解$m_k$，$m_{k+1}$，$m_{k+2}$以及对应的系数 $a_k$，$b_k$，$c_k$，$d_k$。

### 举例
例如，在 $S''(x_0) = S''(x_n) = 0$ 的自然边界条件下，（等会再给个详细的整理过程）
$$
\begin{pmatrix}
1 & 0 & 0 & 0 & \cdots & 0 \\
h_0 & 2(h_0 + h_1) & h_1 & 0 & \cdots & 0 \\
0 & h_1 & 2(h_1 + h_2) & h_2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & h_{n-2} & 2(h_{n-2} + h_{n-1}) & h_{n-1} \\
0 & 0 & \cdots & 0 & 0 & 1
\end{pmatrix}
\begin{pmatrix}
m_0 \\
m_1 \\
m_2 \\
\vdots \\
m_{n-1} \\
m_n
\end{pmatrix}
= 6
\begin{pmatrix}
\frac{y_2 - y_1}{h_1} - \frac{y_1 - y_0}{h_0} \\
\frac{y_3 - y_2}{h_2} - \frac{y_2 - y_1}{h_1} \\
\vdots \\
\frac{y_n - y_{n-1}}{h_{n-1}} - \frac{y_{n-1} - y_{n-2}}{h_{n-2}} \\
0
\end{pmatrix}
$$
可以利用Gauss消元法 / LU分解法 / Gauss-Seidel迭代法求解上述矩阵。

### 代码实现
```python
def 
```

# 数据拟合
数据拟合是通过数学模型近似描述离散数据点的方法，与插值的区别在于：
- 插值：强制通过所有数据点
- 拟合：允许存在偏差，追求整体趋势匹配

基本思想：偏差最小化。

设逼近函数$\phi(x)$和真实数据$y_i$ ($i=1, 2, ..., n$)。则最小化偏差可以由以下方式定义：
$$
L_1\text{范数（绝对偏差和）}：\min \sum_{i=1}^n |\phi(x_i) - y_i|\\
L_2\text{范数（平方偏差和）}：\min \sum_{i=1}^n (\phi(x_i) - y_i)^2\\
L_{\infty}\text{范数（最大偏差）}：\min {\max_{1 \leq i \leq n}|\phi(x_i) - y_i|}
$$

## 最小二乘法拟合
数学原理：**最小化残差平方和**
$$
\min \sum_{i=1}^n (y_i - f(x_i))^2
$$

### 线性拟合
考虑$f(x) = a + bx$。

$\therefore$ 残差平方和：$Q(a,b) = \sum_{i=1}^n (f(x_i) - y_i)^2 = \sum_{i=1}^n (a + b x_i - y_i)^2$

为使 $Q(a, b)$ 最小，需满足：
$$
\begin{cases}
\frac{\partial Q}{\partial a} = 2\sum_{i=1}^n (a + b x_i - y_i) = 0\\
\frac{\partial Q}{\partial b} = 2\sum_{i=1}^n (a + b x_i - y_i)x_i = 0
\end{cases}
$$
将方程组整理为矩阵形式：
$$
\begin{pmatrix}
n & \sum x_i \\
\sum x_i & \sum x_i^2
\end{pmatrix}
\begin{pmatrix}
a \\
b
\end{pmatrix}
=
\begin{pmatrix}
\sum y_i \\
\sum x_i y_i
\end{pmatrix}
$$
求解$(a, b)$即可。（同前，可以用Gauss消元法 / LU分解法 / Gauss-Seidel迭代法）
### 多项式拟合
k次多项式模型：
$$
f(x) = \sum_{j=0}^k a_j x^j = a_0 + a_1 x + a_2 x^2 + ... + a_k x^k
$$
残差表达式
$$
Q(a_0, a_1, ..., a_k) = \sum_{i=1}^n(a_0 + a_1 x_i + a_2 x_i^2 + ... + a_k x_i^k - y_i)^2
$$
对系数求导：
$$
\begin{cases}
\frac{\partial Q}{\partial a_0} = 2 \sum_{i=1}^n(a_0 + a_1 x_i + a_2 x_i^2 + ... + a_k x_i^k - y_i) = 0\\
\frac{\partial Q}{\partial a_1} = 2 \sum_{i=1}^n (a_0 + a_1 x_i + a_2 x_i^2 + ... + a_k x_i^k - y_i) \cdot x_i = 0\\
\frac{\partial Q}{\partial a_2} = 2 \sum_{i=1}^n (a_0 + a_1 x_i + a_2 x_i^2 + ... + a_k x_i^k - y_i) \cdot x_i^2 = 0\\
\vdots\\
\frac{\partial Q}{\partial a_k} = 2 \sum_{i=1}^n (a_0 + a_1 x_i + a_2 x_i^2 + ... + a_k x_i^k - y_i) \cdot x_i^k = 0
\end{cases}
$$
整理：
$$
\begin{cases}
a_0 \sum 1 + a_1 \sum x_i + \cdots + a_k \sum x_i^k = \sum y_i \\
a_0 \sum x_i + a_1 \sum x_i^2 + \cdots + a_k \sum x_i^{k+1} = \sum x_i y_i \\
\vdots \\
a_0 \sum x_i^k + a_1 \sum x_i^{k+1} + \cdots + a_k \sum x_i^{2k} = \sum x_i^k y_i
\end{cases}
$$
整理成矩阵形式：
$$
\begin{pmatrix}
n & \sum x_i & \cdots & \sum x_i^k\\
\sum x_i & \sum x_i^2 & \cdots & \sum x_i^{k+1} \\
\vdots & \vdots & \ddots & \vdots \\
\sum x_i^k & \sum x_i^{k+1} & \cdots & \sum x_i^{2k}
\end{pmatrix}
\begin{pmatrix}
a_0 \\
a_1 \\
\vdots \\
a_k
\end{pmatrix}
=
\begin{pmatrix}
\sum y_i \\
\sum x_i y_i \\
\vdots \\
\sum x_i^k y_i
\end{pmatrix}
$$
求解即可。
### 非线性拟合
| 模型类型 | 变换方法 | 线性形式 |
|---------|----------|----------|
| 指数模型 $y=ae^{bx}$ | 取对数 | $\ln y = \ln a + bx$ |
| 幂函数 $y=ax^b$ | 取对数 | $\ln y = \ln a + b\ln x$ |
| ？对数模型 $y=a\ln x + b$ | - | 直接拟合 |
|倒数模型 $y = \frac{x}{a+bx}$| 取倒数 | $\frac{1}{y} = \frac{a}{x} + b$ |
- 核心思想：参数线性化

容易出现的问题：欠拟合&过拟合。（有机器学习基础的应该都了解……不解释了。）


# 数值微分
利用离散的点进行数值微分。

## Case 1：非边界点微分值的计算
假设数据点均匀分布，$x_i = x_0 + i \cdot h$，在$x_{i-1}$、$x_i$、$x_{i+1}$附近Taylor展开：

$$
f_{i+1} = f_i + h f'_i + \frac{h^2}{2} f''_i + \frac{h^3}{6} f'''_i + O(h^4)
$$
$$
f_{i-1} = f_i - h f'_i + \frac{h^2}{2} f''_i - \frac{h^3}{6} f'''_i + O(h^4)
$$
就可以得到以下公式：

| 类型 | 公式 | 误差阶 |
|------|------|--------|
| 前向差分 | $f'(x) \approx \frac{f_{i+1}-f_i}{h}$ | $O(h)$ |
| 后向差分 | $f'(x) \approx \frac{f_i-f_{i-1}}{h}$ | $O(h)$ |
| 中心差分 | $f'(x) \approx \frac{f_{i+1}-f_{i-1}}{2h}$ | $O(h^2)$ 

再在$x_{i-2}$、$x_{i+2}$处展开：
$$
f_{i+2} = f_i + 2h \cdot f_i' + \frac{4h^2}{2!} f_i'' + \frac{8h^3}{3!} f_i'''+...
$$
$$
f_{i-2} = f_i - 2h \cdot f_i' + \frac{4h^2}{2!} f_i'' - \frac{8h^3}{3!} f_i'''+...
$$
$f_{i-2}+8f_{i+1}-f_{i+2}-8f_{i-1}$ 消去 $f_i'''$、$f_i''$ 和 $f_i$，可得：
$$
f'_i = \frac{f_{i-2}-8f_{i-1}+8f_{i+1}-f_{i+2}}{12h} + O(h^4)
$$
类似地，可以得到三点或五点的二阶微分：
$$
f_i'' = \frac{f_{i+1} - 2f_i + f_{i-1}}{h^2} + O(h^2)
$$
$$
f_i'' = \frac{-f_{i-2} + 16f_{i-1} - 30f_i + 16f_{i+1} - f_{i+2}}{12h^2} + O(h^4)
$$

## Case 2：边界点微分值的计算
以$x_0$为例。

- 利用两点计算一阶微分：
   $$
   f'(x_0) = \frac{f_1-f_0}{h} + O(h)
   $$
   缺陷：精度只有$O(h)$。

- 利用三点计算一阶微分：
   $$
   \begin{cases}
   f_0 = f_0\\
   f_1 = f_0 + h \cdot f_0' + \frac{h^2}{2!} \cdot f_0'' + O(h^3)\\
   f_2 = f_0 + 2h \cdot f_0' + \frac{4h^2}{2!} \cdot f_0'' + O(h^3)
   \end{cases}
   $$
   通过他们的线性组合将 $f_0''$ 消去，我们可以有：
   $$
   f_0' = \frac{4f_1 - 3f_0 - f_2}{2h} + O(h^2)
   $$

同理，$x_n$处利用三点计算一阶微分：
$$
f_n' = \frac{3f_n + f_{n-2} - 4f_{n-1}}{2h} + O(h^2)
$$

# 数值积分
实际中常见情况：$f(x)$ 表达式未知，仅以离散点列 $(x_i,f(x_i))$ 给出 or 很难求出 $f(x)$ 的原函数。

$\Rightarrow$ 可以用分割小区间再数值求和的方法求解积分。

由前面知识可知，拉格朗日插值$L_n(x) = \sum_{k=0}^n l_k(x) \cdot f(x_k)$，其中 $l_k(x)$为插值基函数。用 $L_n(x)$ 代替 $f(x)$ 可得：
$$
I(f) = \int_a^b f(x) \text{d}x \overset{f(x)换成L_n(x)}{=} \int_{a}^{b} L_n(x) \text{d}x = \int_{a}^{b} \sum_{k=0}^{n} l_k(x) f(x_k) \text{d}x
$$
这里的积分与求和可以交换顺序（别问我为什么……我不是学数学的）：
$$
\int_{a}^{b} \sum_{k=0}^{n} l_k(x) f(x_k) \text{d}x = \sum_{i=0}^{n} \left[ \int_{a}^{b} l_i(x) dx \right] f(x_i) = \sum_{i=0}^{n} \alpha_i f(x_i)
$$
此时 $x_i$ 是求积分节点，$\alpha_i$ 是求积分系数。

**Newton-Cotes积分方法：** 对积分区间 $[a,b]$ 作 $n$ 等分，步长 $h = \frac{b-a}{n}$ 进行积分节点构造。

## 非复化数值积分
下面就以Newton-Cotes积分方法为例，查看一阶、二阶和多阶插值下的计算情况。
1. 一阶梯形积分：我们知道两点插值 ($n=1$)
   $$
   L_1(x) = l_0(x)y_0 + l_1(x)y_1 = \frac{x-x_1}{x_0-x_1} \cdot y_0 + \frac{x-x_0}{x_1-x_0} \cdot y_1
   $$
   代回积分表达式：
   $$
   \begin{aligned}
   \int_a^b f(x) \text{d}x = \int_a^b L_1(x) \text{d}x &= \int_a^b \sum_{i=0}^1 l_i(x) f(x_i) \text{d}x\\
   &= \int_a^b l_0(x) \text{d}x \cdot f(a) + \int_a^b l_1(x) \text{d}x \cdot f(b)\\
   &= \frac{b - a}{2} \left[ f(a) + f(b) \right]
   \end{aligned}
   $$
   【 这里计算一下$\int_a^b l_0(x) \text{d}x$：
   $$
   \int_a^b l_0(x) \text{d}x = \int_a^b \frac{x-b}{a-b} \text{d}x = \frac{1}{a-b} (\frac{x^2}{2}-bx)\mid^b_a = \frac{1}{2(a-b)}(b-a)(a-b)
   $$
   $\int_a^b l_1(x) \text{d}x$同理 】

   误差公式前面也已得出：
   $$
   R(x) = f(x) - P_1(x) = \frac{f''(\xi)}{2}(x-x_0)(x-x_1)
   $$
   代入积分：
   $$
   \begin{aligned}
   E_1(x) &= \int_a^b \frac{f''(\xi)}{2}(x - a)(x - b) \text{d}x\\
   &= \frac{f''(\xi)}{2} \int_a^b (x - a)(x - b) \text{d}x = -\frac{f''(\xi)}{12}(b - a)^3
   \end{aligned}
   $$

2. 二阶辛普森（Simpson）积分：$n=2$，三点插值
   $$
   \begin{aligned}
   L_2(x) &= l_0(x) \cdot y_0 + l_1(x) \cdot y_1 + l_2(x) \cdot y_2\\
   &= \frac{(x - x_1)(x - x_2)}{(x_0 - x_1)(x_0 - x_2)}y_0 + \frac{(x - x_0)(x - x_2)}{(x_1 - x_0)(x_1 - x_2)}y_1 + \frac{(x - x_0)(x - x_1)}{(x_2 - x_0)(x_2 - x_1)}y_2
   \end{aligned}
   $$
   带入积分表达式：
   $$
   \begin{aligned}
   \int_a^b f(x) \text{d}x = \int_a^b L_2(x) \text{d}x &= \int_a^b \sum_{i=0}^2 l_i(x) f(x_i) \text{d}x\\
   &= \int_a^b l_0(x) \text{d}x \cdot f(a) + \int_a^b l_1(x) \text{d}x \cdot f\left(\frac{a+b}{2}\right) + \int_a^b l_2(x) \text{d}x \cdot f(b)\\
   &= \frac{b-a}{6} \left[ f(a) + 4f\left(\frac{a+b}{2}\right) + f(b) \right]
   \end{aligned}
   $$
   （这里的误差项我不确定课件上的是否正确……[ac01]）

   以此类推，可总结如下：（我懒得一一去算了……有空来填这个坑）
   | Step size \( h \)       | Common name         | Formula                                      | Error term                     |
   |-------------------------|---------------------|----------------------------------------------|--------------------------------|
   | $b - a$             | Trapezoidal rule    | $\frac{h}{2}(f_0 + f_1)$                   | $-\frac{1}{12}h^3f''(\xi)$   |
   |  $\frac{b - a}{2}$   | Simpson's rule      | $\frac{h}{3}(f_0 + 4f_1 + f_2)$            | $-\frac{1}{90}h^5f^{(4)}(\xi)$ |
   | $\frac{b - a}{3}$   | Simpson's 3/8 rule  | $\frac{3h}{8}(f_0 + 3f_1 + 3f_2 + f_3)$    | $-\frac{3}{80}h^5f^{(4)}(\xi)$ |
   | $\frac{b - a}{4}$   | Boole's rule        | $\frac{2h}{45}(7f_0 + 32f_1 + 12f_2 + 32f_3 + 7f_4)$ | $-\frac{8}{945}h^7f^{(6)}(\xi)$ |


## 复化数值积分
Runge现象：在高次多项式插值中，随着插值节点数量的增加，插值结果在区间端点附近出现剧烈振荡。（前面已经提过）

$\therefore$ 实际计算过程（复化数值积分思想）：把积分区间分割成多个子区间，然后在每个子区间采用一阶梯形积分或二阶辛普森积分，再将所有子区间积分值加起来。

1. 复化梯形积分：
   对被积函数 $f(x)$ 在积分区间 $[a, b]$ 内作等距分割，$h = \frac{b-a}{n}$，并且在每个区间内做线性插值，那么
   $$
   \begin{aligned}
   \int_{a}^{b} f(x) \text{d}x = \sum_{i=0}^{n} \int_{x_i}^{x_{i+1}} f(x) \text{d}x &= \sum_{i=0}^{n} \int_{x_i}^{x_{i+1}} L_1(x) + R_n(x) \text{d}x\\
   &= \sum_{i=0}^{n} \left( \frac{h}{2} [f(x_i) + f(x_{i+1})] - f''(\eta_i) \frac{h^3}{12} \right)\\
   &= \sum_{i=0}^{n} \left( \frac{h}{2} [f(x_i) + f(x_{i+1})] - f''(\eta_i) \frac{h^3}{12} \right)\\
   &\approx h \left[ \frac{1}{2} f(a) + \sum_{i=1}^{n-1} f(x_i) + \frac{1}{2} f(b) \right]
   \end{aligned}
   $$
   同时我们可以估计其误差为（第二个等号到第三个等号，以及最后两个等号之间均用了积分中值定理）：
   $$
   \begin{aligned}
   E_n(f) = \sum_{i=0}^{n} \left( -f''(\eta_i) \frac{h^3}{12} \right) &= \frac{h^2}{12}\sum_{i=0}^{n} (-f''(\eta_i) \cdot h) \\
   &= -\frac{h^2}{12}\sum_{i=0}^{n}\int_{x_i}^{x_{i+1}} f''(x) \text{d}x \\
   &= -\frac{h^2}{12} \int_{a}^{b} f''(x) \text{d}x\\
   &= -\frac{h^2}{12} (b-a) f''(\eta)\\
   \end{aligned}
   $$
   又 $\because$ $h = \frac{b-a}{n} \;$ 代入上式：
   $$
   \text{L.H.S.}= -\frac{(b-a)^3}{12n^2} f''(\eta), \quad \eta \in (a, b)
   $$

2. 复化 Simpson 积分：
   对各子区间做二阶多项式插值。

   对积分区间 $[a,b]$ 内作偶数等距分割，$n=2m$，此时我们有 $2m+1$ 个积分节点和 $m$ 个积分区间，积分步长 $h = \frac{b-a}{n}$。每次积分在 $(x_{2i}, x_{2i+1}, x_{2i+2})$三点上进行：
   $$
   \begin{aligned}
   \int_a^b f(x) \text{d}x &= \sum_{i=0}^{m-1} \int_{x_{2i}}^{x_{2i+2}} f(x) \text{d}x\\
   &= \sum_{i=0}^{m-1} \{ \frac{2h}{6} \left[ f(x_{2i}) + 4f(x_{2i+1}) + f(x_{2i+2}) \right] - \frac{(2h)^5}{2880} f^{(4)}(\eta_i) \} (\text{这里步长为}2h)\\
   &= \frac{h}{3} \left[ f(a) + 4 \sum_{i=0}^{m-1} f(x_{2i+1}) + 2 \sum_{i=1}^{m-1} f(x_{2i}) + f(b) \right] + E_n(f)
   \end{aligned}
   $$
   其中误差项 
   $$
   \begin{aligned}
   E_n(f)=\sum_{i=0}^{m-1}-\frac{(2h)^5}{2880} f^{(4)}(\eta_i) &\approx -\frac{h^4}{180}(b-a)f^{(4)}(\eta)\\
   &= -\frac{(b-a)^5}{180n^4}f^{(4)}(\eta) \quad \eta \in (a,b)
   \end{aligned}
   $$
   （约化技巧与复化梯形积分误差项类似，注意这里的步长为$2h$。）

# 数据滤波
实际观测中永远免不了随机噪声。噪声分类：
- 过程噪声：物理系统演化过程中系统本身的随机性
- 测量噪声：做测量时测量行为带来的随机性
## 滤波方法一：平滑滤波
最简单的平滑滤波：移动窗口的平均法。

对于一组等间隔时间序列数据：$y_0, y_1, \cdots, y_n$，我们可以选择一个数据窗口$\pm k$，对窗口内数据进行平均化：
$$
\hat y_i = \frac{1}{2k+1} \sum_{-k}^{k} y_{i+k}
$$
统计涨落正比于$\frac{1}{\sqrt{2k+1}}$。

{% note info simple %}

统计涨落的计算:

假设噪声数据 $\{\epsilon_i\}$ 独立同分布，均值为0，方差为$\sigma$。

对窗口内$2k+1$个点取平均：$\hat y_i = \frac{1}{2k+1} \sum_{-k}^{k} y_{i+k}$，则噪声的方差等于：
$$
\text{Var}(\hat \epsilon_i) = \text{Var}(\frac{1}{2k+1} \sum^{k}_{-k}\epsilon_{i+k})
$$

根据方差计算性质，后者 $= \frac{1}{(2k+1)^2} \text{Var}(\sum^{k}_{-k}\epsilon_{i+k}) = \frac{\sigma}{(2k+1)^2} \cdot (2k+1) = \frac{\sigma}{2k+1} $。

因此，标准差（即统计涨落）正比于$\frac{1}{\sqrt{2k+1}}$。

{% endnote %}

加权平滑滤波（Savitzky-Golay Filter）：在一个移动窗口内，中心处的权重应该更大，这样可以更好地保留信息。



## 滤波方法二：卡尔曼滤波