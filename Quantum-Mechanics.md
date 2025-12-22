
# Quantum-Mechanics
量子力学について，プログラミングを交えて学ぶリポジトリ

## シュレーディンガー方程式の離散化
1. 復習
    1. 粒子の位置や運動量は確定値を持たない->波動関数 $\psi (x)$ によって確率分布として記述される．
    2. 物理量は演算子として作用->ある物理量Aを観測するとは，その演算子 $\hat A$ の固有値問題を解くことと同義．
  - $\hat A \psi = a \psi$ が成り立つときのみ，確定した値aが観測できる．

2. シュレーディンガー方程式は時間に依存しない
- 1次元のポテンシャル $V(x)$ 中の粒子を考えると，エネルギー固有状態を求める方程式は以下の通り．

```math
-\frac{\hbar ^2}{2m}\frac{d^2\psi (x)}{dx^2}+V(x)\psi (x) = E\psi (x)
```
- どうやって計算機的に解こう？
- 実は上の方程式は $H\psi = E\psi$ のような固有値問題
  - 左辺の微分演算子とポテンシャル項を合わせてハミルトニアン $\hat H$

3. 空間の離散化(差分法)
- 空間xを刻み幅 $\Delta x$ で離散化-> $\psi (x)$ は配列になる．
- 二階微分を差分近似すると？

```math
\frac{d^2 \psi}{dx^2} \approx \frac{\psi_{i+1}-2\psi_{i}+\psi_{i-1}}{(\Delta x)^2}
```
- これを方程式に代入して，定数項 $\frac{\hbar^2}{2m}$ を省略して整理すると漸化式になる．

```math
-\frac{\psi_{i-1}}{(\Delta x)^2}+(\frac{2}{(\Delta x)^2}+V_{i})\psi_{i}-frac{\psi{i-1}}{(\Delta x)^2}=E\psi_{i}
```
- 行列形式に直せば
```math
\begin{pmatrix} \ddots & \ddots & 0 \\ -1 & 2+v_{i}(\Delta x)^2 & -1 \\ 0 & \ddots & \ddots \end{pmatrix} \begin{pmatrix} \vdots \\ \psi_{i} \\ \vdots \end{pmatrix} = E(\Delta x)^2 \begin{pmatrix} \vdots \\ \psi_{i} \\ \vdots \end{pmatrix}
```

- ベクトル $\psi$ に対して行列 $H$ がかかっている．
- 上の**三重対角行列**をプログラムに実装したい->Create-Hamiltonian
- 条件
  - $\psi (x)->0\le x\le N$
  - $V(x)->xの範囲内で0, それ以外で\infty$
  - ディリクレ境界条件-> $\psi (0)=0$ , $\psi (N)=0$

- ディリクレ境界条件はなぜ必要？
- 方程式:
```math
-\psi_{i-1} + (2 + V)\psi_{i} - \psi_{i+1} = E\psi_{i}
```
- i=1,Nのときにこの境界条件が発動する．
- 本来 $1\le i\le N$ だから，方程式にi=1, Nを代入したときに現れる $\psi_{0}, \psi_{N+1}$ は未知数．
- ディリクレ境界条件はこの未知数を前もって定めておくためのもの．
  - これによって，計算は範囲内のものだけで良い，ということになる．
### 行列の固有値を求めるアルゴリズム
#### 固有値の復習
- 定義:
```math
\mathbf{A}\vec{x}=\lambda \vec{x}
```
- このとき $\vec{x}$ を $\mathbf{A}$ の固有ベクトル(eigenvector)， $\lambda$ を固有値(eigenvalue)という．
  - $\mathbf{A}$ は正方行列．固有ベクトルは $\vec{0}$ でない． $\lambda$ はスカラ．

- 固有値すべてを集めた集合を**固有空間**という．
- 固有値の求め方
  - $det(\mathbf{A}-\lambda I)=0)$ を解けばいい．
  - 何故か？
    - とりあえず行列式とはなにかから考えてみよう．

- 行列式の定義: $\mathbf{A}$ が $N\times N$ 行列のとき
```math
\begin{equation}
det\mathbf{A} = \sum_{\sigma \in S_{n}}sgn(\sigma)\prod_{i=1}^{n}a_{i\sigma (i)}
\end{equation}
```
- $3\times 3$ のときで確認してみる．
1. $\mathbf{A}$ を
```math
\mathbf{A}=\begin{pmatrix} a_{11} & a_{12} & a_{13}\\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a{33} \end{pmatrix}
```
とする.

2. $det\mathbf{A}$ は以下で表される．ただし，組み合わせは $3!=6$ 通りあるので, $\sigma_{n}$ は６つ．
```math
\begin{matrix}
det\mathbf{A}=\sum_{\sigma \in S_{3}}sgn(\sigma)\prod_{i=1}^{3}a_{i\sigma (i)} \\
=sgn(\sigma_{1})a_{1\sigma_{1}(1)}a_{2\sigma_{1}(2)}a_{3\sigma_{1}(3)} \\
+sgn(\sigma_{2})a_{1\sigma_{2}(1)}a_{2\sigma_{2}(2)}a_{3\sigma_{2}(3)} \\
\vdots \\
+sgn(\sigma_{6})a_{1\sigma_{6}(1)}a_{2\sigma_{6}(2)}a_{3\sigma_{6}(3)}
\end{matrix}
```

3. sgn記号は置換の回数で符号を決める．奇置換->-1, 偶置換->1
```math
\sigma_{1}=
\begin{bmatrix}
1 & 2 & 3 \\
1(ここが\sigma_{1}(1)) & 2(\sigma_{1}(2)) & 3(\sigma_{1}(3))
\end{bmatrix} 
```
- よって，1項目は $a_{11}a_{22}a_{33}$  .偶置換であることに注意
- 一度置換して
```math
\sigma_{2}=
\begin{bmatrix} 1&2&3 \\ 2&1&3 \end{bmatrix}
```
- 2項目は $-a_{12}a_{21}a_{33}$ .奇置換であることに注意．
- 同じように繰り返していく
```math
\sigma_{3}=
\begin{bmatrix} 1&2&3 \\ 3&1&2 \end{bmatrix}
```
```math
\sigma_{4}=
\begin{bmatrix} 1&2&3 \\ 3&2&1 \end{bmatrix}
```
```math
\sigma_{5}=
\begin{bmatrix} 1&2&3 \\ 2&3&1 \end{bmatrix}
```
```math
\sigma_{6}=
\begin{bmatrix} 1&2&3 \\ 1&3&2 \end{bmatrix}
```
- 結果はサラスの公式を使ったときと同じものになる．

- 重要な性質
    - $N\times N$ 行列Aをn個の列ベクトルで表現した
```math
\mathbf{A}=\begin{pmatrix} a_{1} & a_{2} & ... & a_{n} \end{pmatrix}
```
には以下のことが成り立つ．
1. 各列ベクトルが線形独立である場合 $det\mathbf{A}\neq 1$ である．
2. 各列ベクトルが線形従属である場合 $det\mathbf{A}=0$ である．->これが答えになる．
- 線形独立とは
    - ベクトル $\mathbf{x_{n}}$ に対して， 
```math
c_{1}\mathbf{x_{1}}+...+c_{n}\mathbf{x_{n}} = 0
```
を満たす係数 $c_{i}$ がすべて0であるとき，ベクトル $x_{n}$ を線形独立という．<br>

- 線形従属とは
    - ベクトル $\mathbf{x_{n}}$ に対して， 
```math
d_{1}\mathbf{x{1}}+...+d_{n}\mathbf{x_{n}}=0
```
を満たす係数 $d_{i}$ に0で無いものが存在するならば，ベクトル $x_{n}$ を線形従属という．

- $(A-\lambda I)\mathbf{x}=0$ に置いて $\mathbf{x}\neq 0$ ならば， $(A-\lambda I)=0$ となるしかなく，このときその列ベクトルたちは線形従属(一次従属)となる．

#### べき乗法(Power Iteration)
- n次正方行列の最大固有値と最大固有値に対応する固有ベクトルを計算する．


#### 逆べき乗法(Inverse Power Iteration)
