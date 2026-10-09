# SQP 的二阶结构：Lagrangian Hessian、Reduced Hessian 与步分解

> **核心结论**：约束最优点的曲率应沿可行流形观察。Lagrangian Hessian 把约束弯曲的影响计入；reduced Hessian 把这个曲率限制到允许运动的切空间。

本文接续第一篇，采用 $L=f+\lambda^Tc+\mu^Tg$。以下几何推导主要针对正则等式约束；不等式的二阶条件最后单独说明。

## 1. $B_k$ 的来源与更新次序

对 $x$ 求二阶导时固定乘子：

$$
\nabla^2_{xx}L=\nabla^2f+\sum_i\lambda_i\nabla^2c_i+\sum_j\mu_j\nabla^2g_j.
$$

Exact SQP 取该矩阵，quasi-Newton 用它的近似。乘子把约束曲率的影响加权到优化模型中，但 QP 约束本身仍是线性的。

迭代的因果次序是：先用当前 $(x_k,\lambda_k,\mu_k)$ 构造 $B_k$，再解 QP 得到方向与新乘子估计；接受步骤后得到新点和乘子，下一轮才构造 $B_{k+1}$。KKT 系统不是在这一轮“凭空解出 $B_k$”。

## 2. 为什么不是只用 $\nabla^2f$？

在一个等式 KKT 点 $x_*$ 取可行曲线 $x(t)$，满足 $c(x(t))=0$。设 $v=x'(0)$、$a=x''(0)$。对约束求一次导数：

$$
A_*v=0.
$$

对第 $i$ 个约束再求一次导数：

$$
\nabla c_i^Ta+v^T\nabla^2c_iv=0.
$$

目标沿曲线的二阶导为

$$
\frac{d^2}{dt^2}f(x(t))\bigg|_0=v^T\nabla^2fv+\nabla f^Ta.
$$

第二项是曲线转弯的贡献。由 KKT 的 $\nabla f=-A_*^T\lambda_*$，得到

$$
\nabla f^Ta=\sum_i\lambda_{*,i}v^T\nabla^2c_iv,
$$

因此

$$
\boxed{\frac{d^2}{dt^2}f(x(t))\bigg|_0=v^T\nabla^2_{xx}L_*v.}
$$

这个推导使用了 KKT stationarity；不能对任意点、任意乘子直接宣称二者相等。$\nabla^2f$ 只描述直线方向的曲率，忽略了保持可行时路径必须转弯的影响。

## 3. 三次约束实例：两轮 exact SQP

采用原学习总结中的例子：

$$
f(x,y)=(x-2)^2+(y-1)^2,\qquad c(x,y)=x^3+y-1=0.
$$

各导数为

$$
q=\begin{bmatrix}2(x-2)\\2(y-1)\end{bmatrix},\quad
A=\begin{bmatrix}3x^2&1\end{bmatrix},\quad
\nabla^2c=\begin{bmatrix}6x&0\\0&0\end{bmatrix},
$$

$$
B=\nabla^2_{xx}L=\begin{bmatrix}2+6\lambda x&0\\0&2\end{bmatrix}.
$$

### 3.1 第一轮

从 $(x_0,y_0)=(1,0)$、$\lambda_0=0$ 开始。此时可行，$q_0=(-2,-2)^T$、$A_0=(3,1)$、$B_0=2I$。

$$
\min_{u,v}\ -2u-2v+u^2+v^2,\qquad 3u+v=0.
$$

KKT 方程为

$$
2u+3\widehat\lambda_0=2,\quad 2v+\widehat\lambda_0=2,\quad 3u+v=0.
$$

解得 $d_0=(-0.2,0.6)^T$、$\widehat\lambda_0=0.8$。为演示局部模型，先假定接受完整步和完整乘子更新：

$$
(x_1,y_1)=(0.8,0.6),\qquad \lambda_1=0.8,\qquad c_1=0.112.
$$

完整步破坏了真实约束，尽管它满足线性化约束。实际求解时必须经过第三、四篇的接受测试；这里的完整步不是对所有算法都必然成立的结论。

### 3.2 第二轮

现在

$$
q_1=\begin{bmatrix}-2.4\\-0.8\end{bmatrix},\quad A_1=\begin{bmatrix}1.92&1\end{bmatrix},\quad B_1=\begin{bmatrix}5.84&0\\0&2\end{bmatrix}.
$$

第二轮系统为

$$
\begin{bmatrix}5.84&0&1.92\\0&2&1\\1.92&1&0\end{bmatrix}
\begin{bmatrix}u\\v\\\widehat\lambda_1\end{bmatrix}
=\begin{bmatrix}2.4\\0.8\\-0.112\end{bmatrix}.
$$

数值解为

$$
d_1\approx\begin{bmatrix}0.032840882\\-0.175054493\end{bmatrix},\qquad \widehat\lambda_1\approx1.150108985.
$$

再假定接受完整步，有

$$
(x_2,y_2)\approx(0.832840882,0.424945507),\qquad c_2\approx0.002623876,
$$

$$
B_2\approx\begin{bmatrix}7.747146687&0\\0&2\end{bmatrix}.
$$

这两轮数值已按同一 KKT 系统复算。不要把第 $k$ 轮 QP 的乘子下标与迭代更新后的乘子下标混用。

| 迭代点 | $(x_k,y_k)$ | 当前乘子 | 真实约束残差 |
| --- | --- | --- | --- |
| $k=0$ | $(1,0)$ | $0$ | $0$ |
| $k=1$ | $(0.8,0.6)$ | $0.8$ | $0.112$ |
| $k=2$ | $(0.832840882,0.424945507)$ | $1.150108985$ | $0.002623876$ |

## 4. 用消去约束验证曲率

由 $y=1-x^3$，目标变成

$$
\phi(x)=(x-2)^2+x^6,\qquad \phi''(x)=2+30x^4.
$$

选择切向基底

$$
Z=\begin{bmatrix}1\\-3x^2\end{bmatrix},\qquad AZ=0.
$$

在 KKT 点，$y$ 方向 stationarity 给出 $\lambda=2x^3$，于是

$$
Z^T\nabla^2_{xx}LZ=2+6\lambda x+18x^4=2+30x^4=\phi''(x).
$$

只用 $\nabla^2f=2I$ 会得到 $2+18x^4$，漏掉约束曲率贡献 $12x^4$。

**适用边界**：上式在采用对应乘子关系时成立；任意迭代乘子不一定满足该关系。消元变量 $x$ 也不是弧长，因此曲率的数值依赖参数化，但二阶正性仍有清楚意义。

此外 $\phi''(x)>0$ 对所有实数成立，所以消元后的目标严格凸；这个特例有唯一全局最小点。一般非凸 NLP 没有这一保证。

## 5. Reduced Hessian：只看允许的方向

在等式正则点，$Z\in\mathbb R^{n\times(n-m)}$ 的列张成 $\ker A$，任意切方向为 $v=Zp$。于是

$$
v^TBv=p^T(Z^TBZ)p,\qquad H_R=Z^TBZ.
$$

在等式 KKT 点、LICQ 和适当光滑性下，真实矩阵满足

$$
Z^T\nabla^2_{xx}L_*Z\succ0
$$

是严格局部极小的二阶充分条件。完整 Hessian 在约束禁止的方向上有负曲率，并不自动否定约束极小性。基底可以不正交；改变基底会改变数值表示，但可逆换基不改变正定性。

## 6. Null-space 解法与 normal / tangential step

若当前点可行，令 $d=Zp$，等式 QP 简化为

$$
(Z^TBZ)p=-Z^Tq.
$$

若当前点不可行，先取 $d_N$ 满足 $Ad_N=-c$，再写 $d=d_N+Zp$：

$$
(Z^TBZ)p=-Z^T(q+Bd_N).
$$

这与完整 KKT 系统在相应非奇异条件下等价。求解时用线性求解、QR 或 SVD，不需要显式求逆。

满行秩时，最小范数 normal step 可以写成

$$
d_N=-A^T(AA^T)^{-1}c,
$$

其职责是修复一阶可行性；tangential step $d_T=Zp$ 满足 $Ad_T=0$，利用剩余自由度改善优化模型。

远离解时，normal step 可能只能部分修复残差，例如在信赖域里解 $\min\|c+Ad_N\|_2^2/2$。切向步还要遵守总步长预算。该分解是几何和算法工具，不要求所有 SQP 实现都显式分两次求解。

## 7. KKT 残差与二阶判别不是一回事

考虑 $\min x$，约束 $x^2+y^2=1$。在 $(1,0)$ 取 $\lambda=-1/2$，有 $c=0$、$\nabla_xL=0$，却是最大点。取 $Z=(0,1)^T$，有

$$
Z^T\nabla^2LZ=-1.
$$

在 $(-1,0)$，乘子为 $1/2$，reduced Hessian 为 $1$，才是极小点。

> **核心结论**：可行残差小说明近似可行，stationarity 残差小说明近似一阶驻定；真正的切向曲率用于区分局部性质。即使满足局部极小条件，也不保证一般非凸问题的全局最优。

对不等式，二阶条件应在**临界锥**上检查；弱活跃约束不能一律当作等式。只有在适当的活跃集与严格互补等条件下，才可简化为活跃约束 Jacobian 的 null-space 检查。

## 8. 第二篇速查

$$
B\approx\nabla^2_{xx}L,\qquad AZ=0,\qquad H_R=Z^TBZ,\qquad d=d_N+Zp.
$$

**常见误区**：把 $B$ 当成只含目标的 Hessian；把近似 BFGS 的正定性当成真实二阶证明；在不可行点只走 $Zp$；把 stationarity 当成全局最优；忽略曲线推导里的 KKT 前提。

下一篇讨论：局部模型给了一个方向后，如何决定接受、缩步、修正或恢复可行性。
