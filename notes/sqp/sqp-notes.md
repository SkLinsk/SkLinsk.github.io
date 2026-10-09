# SQP、Filter 与 Funnel：四篇学习笔记

按照 SQP 基础 → 二阶结构 → Filter → Funnel 的顺序阅读。

---

# SQP 基础：从约束线性化到 KKT–Newton

> **核心结论**：SQP 每轮用一个二次规划产生候选方向。它同时处理目标下降与约束的一阶修复；线性化可行并不意味着真实非线性可行。

本文是四篇系列的第一篇。先理解方向如何产生，再进入第二篇的二阶结构、第三篇的 Filter 和第四篇的 Funnel。公式统一采用列向量梯度、按行排列的约束 Jacobian，以及 $g(x)\le0$ 的不等式约定。

## 1. SQP 解决什么问题？

考虑光滑非线性规划（NLP）：

$$
\begin{aligned}
\min_{x\in\mathbb R^n}\quad &f(x),\\
\text{s.t.}\quad &c(x)=0,\quad g(x)\le0.
\end{aligned}
$$

其中 $c:\mathbb R^n\to\mathbb R^m$，$g:\mathbb R^n\to\mathbb R^p$。负梯度只告诉我们如何降低目标，不知道哪些方向会破坏约束；可行域通常还是弯曲的。

SQP 的策略是：在当前点把约束线性化，把优化结构近似为二次模型，求出方向后再决定走多远。

$$
\text{NLP}\longrightarrow \text{QP}_0\longrightarrow \text{QP}_1\longrightarrow\cdots,
\qquad x_{k+1}=x_k+\alpha_kd_k.
$$

“序列”意味着在新点重建模型，不是反复求解一个固定 QP。

## 2. Jacobian：变量移动如何改变约束？

把各约束梯度转置后堆成行：

$$
A_k=J_c(x_k)=\begin{bmatrix}\nabla c_1(x_k)^T\\\vdots\\\nabla c_m(x_k)^T\end{bmatrix}\in\mathbb R^{m\times n}.
$$

一阶展开为

$$
c(x_k+d)=c_k+A_kd+O(\|d\|^2),\qquad c_k=c(x_k).
$$

$c_k$ 是当前残差，$A_kd$ 是沿方向移动时约束的一阶变化。于是 $c_k+A_kd=0$ 表示：**局部线性模型预测下一点可行**。

若 $c_k=0$ 且 $A_k$ 满行秩，$A_kd=0$ 描述可行流形的切空间；若 $c_k\ne0$，则 $A_kd=-c_k$ 是修复残差的仿射条件。

> **容易混淆的地方**：可行点的切向条件是 $Ad=0$，不可行点的 SQP 约束通常是 $Ad=-c$。后者不是经过原点的向量空间。

## 3. SQP 的 QP 子问题

记 $q_k=\nabla f(x_k)$、$C_k=J_g(x_k)$，避免把目标梯度与不等式函数 $g$ 混用。

$$
\begin{aligned}
\min_d\quad &q_k^Td+\tfrac12d^TB_kd,\\
\text{s.t.}\quad &c_k+A_kd=0,\\
&g_k+C_kd\le0.
\end{aligned}
$$

常数 $f_k$ 不影响最小解，故省略。理论上的 exact SQP 取

$$
B_k=\nabla^2_{xx}L(x_k,\lambda_k,\mu_k),\quad
L=f+\lambda^Tc+\mu^Tg,\quad \mu\ge0.
$$

它不是只对目标做 Taylor 展开后得到的普通二次模型：$B_k$ 含约束曲率，第二篇会推导其来源。QP 的约束仍然是线性的。

## 4. 单位圆上的第一步

考虑

$$
f(x,y)=(x-2)^2+(y-1)^2,\qquad c(x,y)=x^2+y^2-1=0.
$$

几何上是在单位圆寻找最靠近 $(2,1)$ 的点，真实解为

$$
(x_*,y_*)=\left(\frac2{\sqrt5},\frac1{\sqrt5}\right).
$$

从 $(1,0)$ 出发，取初始乘子 $\lambda_0=0$，于是 $B_0=2I$。当前信息为

$$
q_0=\begin{bmatrix}-2\\-2\end{bmatrix},\qquad
A_0=\begin{bmatrix}2&0\end{bmatrix},\qquad c_0=0.
$$

QP 为

$$
\min_{u,v}\ -2u-2v+u^2+v^2,\qquad 2u=0.
$$

所以 $d_0=(0,1)^T$。它确实沿圆在 $(1,0)$ 的切线移动，但完整步到达 $(1,1)$，真实残差为 $1$。

对任意步长 $\alpha$，有

$$
c(1,\alpha)=\alpha^2,\qquad
f(1,\alpha)=2-2\alpha+\alpha^2.
$$

因此当 $0<\alpha<2$ 时目标下降，约束却出现二阶误差。这不是 QP 算错，而是模型只满足一阶约束。第三篇的接受规则和 SOC 正是为这种差异服务。

## 5. QP 的 KKT 条件

QP 解返回 $d_k$ 和乘子估计 $\widehat\lambda_k,\widehat\mu_k$：

$$
\begin{aligned}
q_k+B_kd_k+A_k^T\widehat\lambda_k+C_k^T\widehat\mu_k&=0,\\
c_k+A_kd_k&=0,\\
g_k+C_kd_k&\le0,\\
\widehat\mu_k&\ge0,\\
\widehat\mu_{k,j}(g_k+C_kd_k)_j&=0.
\end{aligned}
$$

等式情形可以写成 saddle-point 系统：

$$
\begin{bmatrix}B_k&A_k^T\\A_k&0\end{bmatrix}
\begin{bmatrix}d_k\\\widehat\lambda_k\end{bmatrix}
=-\begin{bmatrix}q_k\\c_k\end{bmatrix}.
$$

即使 $B_k\succ0$，整个 KKT 块矩阵通常也是不定矩阵。若 QP 凸，KKT 在适当条件下刻画最优解；若 QP 非凸，仅满足 KKT 不能保证取到最小点。

## 6. SQP 为什么等价于 KKT–Newton？

只考虑等式。原 NLP 的一阶方程为

$$
F(x,\lambda)=\begin{bmatrix}\nabla_xL(x,\lambda)\\c(x)\end{bmatrix}=0.
$$

Newton 线性化得到

$$
\begin{bmatrix}\nabla^2_{xx}L_k&A_k^T\\A_k&0\end{bmatrix}
\begin{bmatrix}d_k\\\Delta\lambda_k\end{bmatrix}
=-\begin{bmatrix}\nabla_xL_k\\c_k\end{bmatrix}.
$$

令 $\widehat\lambda_k=\lambda_k+\Delta\lambda_k$，第一行变成

$$
q_k+B_kd_k+A_k^T\widehat\lambda_k=0.
$$

这正是等式 QP 的 stationarity。因此 exact SQP 的方向可理解为对原 KKT 方程做 Newton。对不等式，局部固定活跃集后可作类似解释，还须处理活跃集改变和互补性。

> **核心结论**：QP 乘子是更新后的乘子估计；Newton 系统里的 $\Delta\lambda$ 是增量。两者相差当前乘子，不能把两个记号混为一谈。

## 7. Exact 与 quasi-Newton / BFGS

Exact SQP 显式计算 Lagrangian Hessian。Quasi-Newton SQP 用梯度差构造近似，以降低二阶导成本。

定义实际接受的位移 $s_k=x_{k+1}-x_k$，使用同一组新乘子构造

$$
y_k=\nabla_xL(x_{k+1},\lambda_{k+1},\mu_{k+1})
-\nabla_xL(x_k,\lambda_{k+1},\mu_{k+1}).
$$

固定乘子是为了主要捕捉变量移动带来的曲率。BFGS 为

$$
B_{k+1}=B_k-\frac{B_ks_ks_k^TB_k}{s_k^TB_ks_k}
+\frac{y_ky_k^T}{s_k^Ty_k}.
$$

当 $B_k\succ0$ 且 $s_k^Ty_k>0$ 时保持正定。曲率条件失效或分母太小时，应阻尼、跳过或重置；大规模问题可以采用 L-BFGS。正定 BFGS 便于求凸 QP，但不能据此证明真实二阶条件成立。

## 8. 第一篇速查

| 组件 | 回答的问题 |
| --- | --- |
| $J_c$ | 变量移动让约束一阶改变多少？ |
| QP | 如何兼顾局部目标曲率和线性约束？ |
| QP KKT | 如何联立求方向和新的乘子估计？ |
| KKT–Newton | 为什么 QP 是原一阶条件的局部解法？ |
| Globalization | 候选方向可以走多远、能否被接受？ |

**不要误读**：完整步不一定可行；QP 乘子不等于乘子增量；$B_k$ 不只是 $\nabla^2f$；KKT 矩阵不因 $B_k$ 正定就整体正定。

下一篇用可行曲线和三次约束实例解释：为什么真正的二阶结构是 $\nabla^2_{xx}L$，以及如何只在允许的方向上观察它。


---

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


---

# Filter-SQP：接受规则、Switching、SOC 与 Restoration

> **核心结论**：SQP 负责产生局部候选方向；globalization 负责检验这一步是否取得可靠进展。Filter 分开管理目标值与违反度，switching 决定何时必须要求目标充分下降。

第一、二篇解释了 QP 与约束曲率。本篇把局部模型接到真实非线性问题上，最后给出接受流程。下一篇再把同一思路连接到 Funnel。

## 1. 方向算出来，为什么还不能直接走？

QP 中 $c_k+A_kd=0$，真实残差却可能是 $O(\|d\|^2)$；目标的实际变化也未必符合局部模型。远离解时，完整步尤其可能走入模型失真的区域。

这里要区分两种分工：

| 分工 | 代表方法 | 作用 |
| --- | --- | --- |
| 接受策略 | Merit、Filter、Funnel | 定义什么叫足够进展 |
| 步骤控制 | Line search、Trust region | 拒绝后缩步或重算局部步骤 |

线搜索试 $x_k+\alpha d_k$，逐渐缩小 $\alpha$；信赖域限制 $\|d\|\le\Delta$，拒绝后通常缩小半径并重新求子问题。不能简单把“缩信赖域”写成“沿同一个方向缩 $\alpha$”。

“Global convergence”通常指从较广起点范围获得到驻点等的收敛性质，不是寻找全局最小值的保证。

## 2. 先统一违反度

一般约束可取

$$
h(x)=\|c(x)\|_1+\|[g(x)]_+\|_1,\qquad [g]_j^+=\max(g_j,0).
$$

$h\ge0$，且 $h=0$ 等价于原始可行。其他论文可能用二范数、平方残差或记号 $\theta$；比较阈值时必须先确认定义和缩放。

## 3. Merit function：兑换成一个标量

传统 $\ell_1$ merit 为

$$
\Phi_\rho(x)=f(x)+\rho h(x),\qquad \rho>0.
$$

沿一个 merit 下降方向检查

$$
\Phi_\rho(x_k+\alpha d)\le\Phi_\rho(x_k)+\sigma\alpha\Phi_\rho'(x_k;d),\quad 0<\sigma<1.
$$

$h$ 不光滑，所以这里使用单侧方向导数；若导数不为负，就不能把右边当作充分下降要求，应先调整罚参数或方向。

等式 QP 满足 $Ad=-c$ 时，$\|c\|_1$ 的方向导数是 $-\|c\|_1$，所以罚项能抵消某些目标上升。在适当正则条件下，超过解处乘子相应对偶范数的参数可产生局部精确罚性质；但过大又可能过分压制目标优化。

> **直觉**：Merit 把“目标改善”和“可行性改善”兑换成同一种货币；罚参数就是兑换率。

## 4. Filter：同时看两个坐标

Filter 保存历史 pair：

$$
\mathcal F=\{(h_i,f_i)\}.
$$

低 $h$、低 $f$ 都更好。候选点对每个条目至少满足一条带余量的改善条件，例如

$$
\forall(h_i,f_i)\in\mathcal F:\quad
h_t\le(1-\gamma_h)h_i\quad\text{或}\quad f_t\le f_i-\gamma_fh_i,
$$

其中 $0<\gamma_h<1$、$\gamma_f>0$。

逻辑是“**对每个历史条目，至少改善一维**”，不是“找到一个条目能改善就算通过”。如果新 pair 两坐标都不大于某旧 pair，则旧 pair 被支配，可删除。

### 4.1 一个可直接判断的数字例子

设 filter 只有 $(h_i,f_i)=(0.5,10)$，$\gamma_h=0.1$、$\gamma_f=0.2$。候选点需满足 $h_t\le0.45$ 或 $f_t\le9.9$。

| 候选 $(h_t,f_t)$ | 对这个条目的判断 | 原因 |
| --- | --- | --- |
| $(0.40,10.4)$ | 通过 | 可行性改善足够，允许目标上升 |
| $(0.60,9.8)$ | 通过 | 目标改善足够，允许违反度上升 |
| $(0.48,9.95)$ | 不通过 | 两项虽然都略有改善，却达不到余量 |
| $(0.55,10.1)$ | 不通过 | 两方面都变差 |

这只是历史 filter 测试，最终接受还需要违反度上限、与当前点的进展检查以及可能触发的目标 Armijo 条件。

## 5. f-type 与 h-type：谁来证明这一步有价值？

**f-type** 是按目标下降验收的步骤；**h-type** 是未进入该分支、按 filter 或可行性规则验收的步骤。

例如 $(f,h):(10,0.5)\to(10.4,0.05)$ 可以通过显著可行性改善获得接受。目标上升不自动代表失败。

但如果 $h\approx10^{-7}$，算法不能长期只依赖极微小的可行性改善而忽略 stationarity。此时需要一个机制要求真正的目标进展。

**严格一点的理解**：某些 filter 算法里的 h-type 是“非 f-type 的可接受步”，不一定每步都严格降低 $h$。类型由算法条件决定，不由名称或事后看到的单个数值决定。

## 6. Switching condition：何时转向目标优化？

教学上可定义二次模型下降

$$
\Delta m_k(d)=-q_k^Td-\tfrac12d^TB_kd,
$$

再比较 $\Delta m_k(d_k)>0$ 且 $\Delta m_k(d_k)\ge\delta h_k^2$。违反度小时阈值小，算法更容易进入 f-type。

不过 $B_k$ 含约束曲率，该量是 QP 模型下降，不是实际目标下降，也不必等于真实 $f$ 的二阶 Taylor 预测。指数 $2$ 更不是所有算法的统一标准。

### 6.1 与线搜索 Armijo 一致的一种规则

令 $D_k=-q_k^Td_k$。本系列完整伪代码采用以下 filter 教学规则：

$$
h_k\le h_{\min},\qquad D_k>0,\qquad
\alpha D_k^{s_f}>\delta h_k^{s_h}.
$$

触发时要求

$$
f_t\le f_k-\sigma\alpha D_k.
$$

因为 $D_k>0$，右边严格低于 $f_k$。如 $f_k=10$、$D_k=2$、$\alpha=1$、$\sigma=0.1$，要求 $f_t\le9.8$；$9.99$ 虽下降，却不充分。

未触发时，仍须对当前点满足 filter 式进展，例如 $h_t\le(1-\gamma_h)h_k$ 或 $f_t\le f_k-\gamma_fh_k$。

这类方向导数规则帮助理解 Wächter–Biegler 的 filter line-search 思路；具体阈值、指数及其他保护条件应按所读算法确定。IPOPT 的算法本体是内点法，不能因为它使用 filter 就整体称为经典 SQP。[算法说明](https://doi.org/10.1007/s10107-004-0559-y)

## 7. Filter 更新：不是接受每个点都插入

本系列约定：h-type 接受后，把**接受前的当前点** $(h_k,f_k)$ 插入 filter，删除被它支配的条目；f-type 不插入。

这样 h-type 的折中成果成为以后必须尊重的历史边界。不同算法可能改变插入点和时机，阅读论文时应把更新规则与接受条件作为一个整体，而不是只复制一条不等式。

## 8. Maratos effect：好的 Newton 步也可能被拒绝

以单条等式为例，若约束有足够光滑性，

$$
c(x_k+d)=c_k+A_kd+\tfrac12d^T\nabla^2c_kd+O(\|d\|^3).
$$

QP 消掉前两项，却留下二阶残差。在解附近，这种残差仍可能使一个本来适合快速局部收敛的完整步被接受策略拒绝，从而反复缩步，妨碍预期的超线性或二次收敛。

这叫 Maratos effect。不是所有完整步拒绝都属于它：方向错误、曲率不可靠或远离解也会导致拒绝。[Filter 局部收敛与 SOC 分析](https://doi.org/10.1137/S1052623403426544)

## 9. SOC：修正二阶约束误差

对完整候选 $x_t=x_k+d_k$，求小修正 $w$：

$$
A_kw\approx-c(x_t),\qquad x_{\mathrm{soc}}=x_k+d_k+w.
$$

在正则、局部适用条件下，$w=O(\|d_k\|^2)$，所以主要修补误差而不替代原 Newton 方向。也可采用更新后的 Jacobian；有不等式时要考虑相关活跃约束并验证其余约束。

单位圆上，从 $(1,0)$ 沿 $d=(0,t)$ 移动，残差是 $t^2$。取 $w=(-t^2/2,0)$ 后

$$
c(1-t^2/2,t)=t^4/4.
$$

它展示误差从二阶降到四阶的机制；真正的 Maratos 判断还要结合目标和接受规则。SOC 候选仍需检查真实函数值，次数与大小要受限，不能未经测试自动接受。

## 10. Restoration：正常优化失效后先恢复可行性

当线性化约束不相容、QP 失败且无法修复，或缩到最小步长仍找不到可接受点时，进入 restoration。

例如优先求解

$$
\min_x\ \tfrac12\|c(x)\|_2^2+\tfrac12\|[g(x)]_+\|_2^2,
$$

也可以采用 elastic/slack 子问题、邻近项或信赖域。暂时把重点从原目标移到违反度，但退出必须满足主算法的重返条件，不能只凭“$h$ 降了一点”就返回。

返回后重新估计乘子，并按规则修复或重置 Hessian 近似。如果停在 $h>0$ 的违反度驻点，应报告恢复失败或不可行驻点；这本身不证明原问题全局无可行解。

Normal step 是正常迭代内部的一阶修复分量；restoration 是正常步骤失效后的独立阶段，两者不能互换名称。

## 11. 接受流程与第三篇速查

```text
QP 方向 → 完整候选
  → 是否满足违反度上限与历史 filter？
  → switching 是否触发？
      f-type：还要目标 Armijo 充分下降
      h-type：还要对当前点取得允许的进展
  → 通过：接受并按类型更新历史
  → 不通过：合适时尝试 SOC，否则缩步
  → 反复失败：restoration，成功后重返主算法
```

**不要误读**：filter 通过不是完整接受测试；h-type 不保证 $h$ 每步单调下降；缩步不自动解决不相容 QP；SOC 不保证可接受；restoration 失败不等于全局不可行；switching 不是逃离局部最优的搜索机制。

下一篇用一个标量上界替代二维历史，解释 Funnel，并把四篇内容汇成完整 SQP 伪代码。


---

# Funnel 与完整 SQP：漏斗边界、恢复机制与复习速查

> **核心结论**：Funnel 维护允许的最大违反度 $\tau_k$，并逐步收紧它。在边界内仍须证明目标或可行性进展；$\tau_k$ 单调不增，不代表 $h_k$ 每步单调下降。

本篇承接 Filter 的步类型与 switching，再给出可逐行复习的 SQP 教学框架。不同文献的接受条件并不完全相同，下面会固定预测量和更新约定，避免混用。

## 1. 从历史边界到漏斗上界

Filter 保存多个 $(h_i,f_i)$；Funnel 主要维护一个违反度上界：

$$
h(x_k)\le\tau_k,\qquad \tau_{k+1}\le\tau_k.
$$

在 $(h,f)$ 平面上，$h>\tau_k$ 的点位于允许区之外。位于边界内只是有资格接受，还要检查相应充分下降条件。

漏斗宽度单调本身不保证 $h\to0$：它可能收敛到正数。完整收敛分析还依赖 switching、目标有界、恢复机制与子问题性质等条件。[统一 funnel restoration SQP 框架](https://arxiv.org/abs/2409.09208)

## 2. 固定一个线搜索教学版本

为与 Armijo 一致，定义 $D_k=-q_k^Td_k$、$P_k(\alpha)=\alpha D_k$。先要求 $h_t\le\tau_k$。

若

$$
D_k>0,\qquad P_k(\alpha)\ge\delta h_k^2,
$$

则按 f-type 验收：

$$
f_t\le f_k-\sigma P_k(\alpha),\qquad \tau_{k+1}=\tau_k.
$$

若不触发该分支，按 h-type 要求

$$
h_t\le\beta\tau_k,\qquad 0<\beta<1,
$$

接受后更新

$$
\tau_{k+1}=(1-\kappa)h_t+\kappa\tau_k,\qquad 0<\kappa<1.
$$

当 $\tau_k>0$ 时，有

$$
h_t\le\tau_{k+1}\le[\kappa+(1-\kappa)\beta]\tau_k<\tau_k.
$$

这既让新点留在更新后的漏斗内，又收窄了上界。某些论文使用二次模型下降 $\Delta m_k$、不同的违反度或其他步类型；替换时必须连同 switching、充分下降与收敛假设一起替换。

## 3. 为什么允许违反度暂时上升？

设 $\tau_0=1$、$h_0=0.7$，一个 h-type 候选为 $h_t=0.4$，取 $\beta=0.8$、$\kappa=0.5$，得到 $\tau_1=0.7$。

下一步若 $h:0.4\to0.5$，仍满足 $h_t<0.7$；只要 f-type 的预测与实际充分下降都通过，就可接受，宽度保持 $0.7$。

| 状态 | 违反度 | 漏斗宽度 | 解释 |
| --- | --- | --- | --- |
| 初始 | $0.7$ | $1.0$ | 位于初始漏斗内 |
| h-type 后 | $0.4$ | $0.7$ | 收窄允许边界 |
| f-type 后 | $0.5$ | $0.7$ | 为目标改善使用边界内的余量 |

这里“h-type”要求相对宽度足够深入，并不严格要求比当前 $h_k$ 更小。Filter 与 Funnel 都不应只凭名称推出逐步单调性。

## 4. Merit / Filter / Funnel 对比

| 方法 | 主要状态 | 对候选点的要求 | 维护任务 |
| --- | --- | --- | --- |
| Merit | 罚参数 $\rho$ | $f+\rho h$ 充分下降 | 处理尺度与罚参数 |
| Filter | 历史 pair 集合 | 对每个条目改善一维，配合当前点测试和 switching | 插入、删除支配条目 |
| Funnel | 违反度上界 $\tau_k$ | 保持在边界内，再满足目标下降或收窄条件 | 保持与更新上界 |

可记成：Merit 设兑换率，Filter 记历史成绩，Funnel 设可行性边界。在某些框架中可用特殊 filter 解释 funnel，但不是任意 funnel 都等于“只存一个历史点的 filter”。

三者都可配合 SOC 与 restoration；线搜索和信赖域则提供拒绝后的步骤控制。

## 5. 停止条件：要检查完整 NLP KKT

采用 $g\le0$、$\mu\ge0$ 的约定，可分别检查

$$
\begin{aligned}
r_s&=\|\nabla_xL(x,\lambda,\mu)\|_\infty,\\
r_e&=\|c(x)\|_\infty,\qquad r_i=\|[g(x)]_+\|_\infty,\\
r_d&=\|[-\mu]_+\|_\infty,\qquad r_c=\|\mu\odot g(x)\|_\infty.
\end{aligned}
$$

各项达到合理缩放后的容差，才报告近似一阶 KKT 收敛。只有等式时，无需不等式对应项。$L$ 的数值本身不是 stationarity residual。

例如 $f(x)=10^{-12}x^2$ 在 $x=10^5$ 时梯度仅为 $2\times10^{-7}$，可能通过宽松绝对容差；因此小残差必须结合尺度解释，不能只看打印出的科学记数法。

步很小或目标变化很小只能作为停滞诊断，不能替代 KKT。小 KKT 残差也不保证局部极小，更不保证一般非凸问题的全局最优。

## 6. 完整 SQP 教学伪代码

这是**模块化线搜索教学框架**，不是某篇论文或生产求解器的逐行实现。Filter 固定采用第三篇的方向导数 switching；Funnel 固定采用本篇的线性预测 $\alpha D_k$。所有恢复退出、曲率修正和 SOC 大小界都须按所选完整算法具体化。

参数：$0<\sigma,r,\beta,\kappa,\gamma_h<1$；$\delta,\gamma_f>0$；步长下界 $\alpha_{\min}>0$。

```text
输入：x0, λ0, μ0≥0, B0, 缩放与容差, 迭代/时间上限
选择接受策略 mode ∈ {MERIT, FILTER, FUNNEL}
初始化 ρ；F=空集；τ0>0 且 τ0≥h(x0)，保留边界余量

for k = 0, 1, ...:
    计算 f, q=∇f, c, g, A=Jc, C=Jg
    若完整 NLP KKT 残差达标：
        返回“一阶 KKT 收敛”及各项残差
    若超过迭代/时间上限：返回“达到上限”

    exact：B=∇²xx L(xk,λk,μk)
    quasi-Newton：使用当前 Bk
    必要时修正曲率或正则化，适配所用 QP 求解策略

    求 QP：min qᵀd + 0.5 dᵀBd
            s.t. c+Ad=0, g+Cd≤0
    得到 d, λhat, μhat
    若 QP 没有有效解：诊断并尝试修复，失败转 RESTORE
    若 d 很小但 KKT 残差仍大：诊断停滞，必要时 RESTORE

    MERIT：调整 ρ，使 φρ'(xk;d)<0
           若无法取得下降方向：RESTORE
    D=-qᵀd；α=1；accepted=false

    while α≥αmin:
        xt=xk+αd，计算真实 ft, ht
        (ok,type)=ACCEPT(xt,α,D,mode)
        若 ok：accepted=true；break

        若 α=1 且拒绝主要源于二阶约束偏差：
            尝试有限次 SOC，求足够小的 w
            xsoc=xk+d+w
            若 xsoc 通过同一真实接受规则：
                xt=xsoc；记录 SOC；accepted=true；break
        α=rα

    若未 accepted：转 RESTORE

    FILTER 且 type=h：将旧点 (hk,fk) 插入 F，删除被支配条目
    FUNNEL 且 type=h：τnew=(1-κ)ht+κτk
    FUNNEL 且 type=f：τnew=τk

    xnew=xt
    选 αdual∈(0,1]，更新乘子估计：
        λnew=λk+αdual(λhat-λk)
        μnew=μk+αdual(μhat-μk)
    或在新点重新估计乘子，保持 μnew≥0

    quasi-Newton：
        s=xnew-xk，包含 SOC 的实际位移
        y=∇xL(xnew,λnew,μnew)-∇xL(xk,λnew,μnew)
        检查曲率与分母，做阻尼 BFGS、跳过或重置
    更新状态并继续

RESTORE:
    解局部违反度最小化或 elastic 子问题
    寻找有足够进展且满足主算法重返条件的 xR
    若成功：
        按算法规则维护 F 或 τ，保证相应不变量
        更新点、重新估计乘子、修复/重置 B，返回主迭代
    否则：报告“恢复失败/不可行驻点/数值失败”等诊断
```

### 6.1 接受模块

```text
ACCEPT(xt,α,D,mode):
    MERIT:
        返回 φρ(xt)≤φρ(xk)+σ α φρ'(xk;d)

    FILTER:
        若 ht 超过允许上限：拒绝
        对 F 中每个 (hi,fi)：
            若 ht≤(1-γh)hi 与 ft≤fi-γf hi 都不成立：拒绝
        switch=(hk≤hmin 且 D>0 且 α D^sf>δ hk^sh)
        若 switch：返回 (ft≤fk-σ αD, f)
        否则返回 (ht≤(1-γh)hk 或 ft≤fk-γf hk, h)

    FUNNEL:
        若 ht>τk：拒绝
        switch=(D>0 且 αD≥δ hk²)
        若 switch：返回 (ft≤fk-σ αD, f)
        否则返回 (ht≤βτk, h)
```

SOC 使用基步的预测量并限制修正大小；具体可证明算法还须明确专用验收细节。Restoration 不能靠随意清空 filter 或放宽 $\tau$ 来绕过进展要求。若 QP 非凸，不能只接受任意驻点方向。

## 7. 复习时最容易混淆的组件

| 组件 | 职责 | 不保证什么 |
| --- | --- | --- |
| SQP / QP | 产生局部候选方向 | 真实非线性可行、完整步可接受 |
| Switching | 选择目标验收或其他进展分支 | 逃离局部极小、全球搜索 |
| SOC | 修补小的二阶约束偏差 | 自动可行或自动通过验收 |
| Restoration | 正常步骤失效时优先降低违反度 | 失败即证明全局无可行解 |
| Reduced Hessian | 检查允许方向的曲率 | 近似矩阵正定就是真实二阶证明 |

## 8. 一页式公式速查

### 模型与方向

$$
L=f+\lambda^Tc+\mu^Tg,\quad \mu\ge0,\quad B\approx\nabla^2_{xx}L.
$$

$$
\min_d q^Td+\tfrac12d^TBd,\qquad c+Ad=0,\quad g+Cd\le0.
$$

### 几何与曲率

$$
Ad_N=-c,\quad AZ=0,\quad d=d_N+Zp,\quad H_R=Z^TBZ.
$$

等式 KKT 点的真实 $H_R\succ0$ 支持严格局部极小；不等式检查临界锥。BFGS 梯度差在新旧两点使用同一组乘子。

### 三种验收与两种补救

**Merit**：$f+\rho h$ 充分下降。**Filter**：对每个历史 pair 至少改善一维，配合 switching / Armijo。**Funnel**：$h_t\le\tau$，f-type 验目标，h-type 收紧上界。

**SOC**：$Aw\approx-c(x+d)$，用小修正消除二阶偏差。**Restoration**：正常优化失效后优先恢复可行性，再满足主算法重返条件。

### 总流程与停止含义

**求导 → 构造 $B$ → 解 QP → 验收 / 缩步 → SOC 或恢复 → 更新点和乘子 → BFGS → 检查完整 KKT。**

小残差是近似一阶条件，真实切向曲率用于局部性质判断；非凸 NLP 一般不保证全局最优。$f$、$h$ 不必每步单调，$\tau$ 单调本身也不保证 $h\to0$。

## 参考与阅读方式

本系列根据提供的 SQP 学习总结重组，保留单位圆与三次约束推导，统一乘子符号，并补充适用条件、数值复核和接受流程。各知识点按“概念 → 数学结构 → 容易卡住的问题”组织。

- [Filter line-search 的局部收敛分析](https://doi.org/10.1137/S1052623403426544)：进一步理解 SOC 与局部收敛。
- [Wächter–Biegler 的实现说明](https://doi.org/10.1007/s10107-004-0559-y)：进一步理解 filter、switching 和 restoration；其算法本体为内点法。
- [A Unified Funnel Restoration SQP Algorithm](https://arxiv.org/abs/2409.09208)：进一步理解漏斗与恢复框架。

不同论文的违反度、预测下降、switching 和历史更新必须一起阅读。本系列伪代码明确标为教学版本，不能直接当作某篇论文的收敛证明对象。
