# Conic Billiard Crypto Lab

有限域上的圓錐曲線反射、撞擊循環與矩陣運算研究玩具。

> **不是可供部署的密碼學函式庫。** 本專案研究「撞擊次數可否作為離散對數形式的秘密」；並未證明任何橢圓／雙曲線實例具有生產等級的密碼學安全性。小質數與示範子群可被直接破解。

## Quick start

Python 3.10+; the core library and tests require no third-party packages.

```bash
python -m unittest discover -s tests -v
python -m conic_billiard scan --p 103 --kind all
python -m conic_billiard.compressed_dh_demo
python tools/verify_large_prime_samples.py
```

The main branch contains mathematical demonstrations and reproducible finite-field
experiments; **it is not approved cryptography**. See [compressed state format](docs/compressed_points.md),
[security notes](docs/security.md) and [large-prime sample results](docs/large_prime_samples.md).

## 研究對象

在奇質數域 \(\mathbb F_p\) 上，統一以

\[
C: F(X)=X^TQX+\ell^TX+c=0
\]

表示圓錐曲線。內建模型：

| 模型 | 方程式 | 地位 |
|---|---|---|
| 橢圓型 | \(x^2/a^2+y^2/b^2=1\) | 主要研究候選 |
| **雙曲線型** | \(x^2/a^2-y^2/b^2=1\) | **已加入共用反射核心、循環搜尋與測試** |
| 正圓型 | \(x^2/r^2+y^2/r^2=1\) | 對照組；可能容易線性化 |
| 拋物線型 | \(y=x^2\) | 對照組；已知特定情形可寫固定 \(2\times2\) 矩陣 |

這裡「橢圓型、雙曲線型」指 **所選方程式的形式**；有限域沒有實數幾何的內外、角度與開放/封閉的完全相同概念。仿射圓錐曲線的等價性也不保證反射映射等價，因為反射還依賴背景雙線性形式。

## 新增：三次橢圓曲線的真正鏡面撞擊（實驗分支）

新增 `conic_billiard/cubic_mirror.py` 與 `docs/cubic_mirror.md`，研究 **Weierstrass 三次曲線**

\[E:y^2=x^3-x\pmod p\]

的鏡面反射，而不是使用既有橢圓曲線點加法。與二次圓錐曲線不同，反射直線通常有 0 或 2 個額外的有限域交點；這是一個**多值代數對應**，尚沒有自然且完整的單值快速倍增或 DH 構造。

```bash
python -m conic_billiard.cubic_experiment --out results/cubic_mirror.json
```

已完整枚舉 `p=101` 的 10,506 個狀態及 `p=1009` 的 1,078,482 個狀態，核對所有有效的反向反射關係，並做了低次多項式不變量搜尋。結果及限制詳見 [三次曲線鏡面反射研究](docs/cubic_mirror.md)。

## 新增：實數封閉分支的真正「第一撞擊點」實驗

前一版的三次曲線測試列舉了反射**直線**的全部交點；它不是實數鏡面撞球的物理落點規則。現在新增 `real_oval.py`，針對 `y²=x³−x` 的實數左側凸形封閉分支，採用**反射射線上的第一個正參數交點**。10,000 組隨機起始狀態全部得到唯一的下一撞擊點，單步反向核對最大誤差約 `1.6e-13`。

```bash
python -m conic_billiard.real_oval_experiment
# 可選：pip install mpmath
python -m conic_billiard.real_oval_high_precision
```

**重要限制：** 實數的「第一個」、「正方向」和「左側封閉連通分支」無法原封不動移植到有限域。這是需要獨立解決的分支選取問題；先前的有限域多值統計仍屬代數交點統計，不能解讀為實數物理撞球的多條後繼。詳見 [第一撞擊點實驗](docs/real_oval_first_hit.md)。

## 對所有圓錐曲線共用的反射公式

對有序撞擊點 \(A,B\in C(\mathbb F_p)\)，令 \(v=B-A\)，並令法向量 \(N=QB+\ell/2\)。如果相關分母可逆：

\[
H_B=I-2\frac{NN^T}{N^TN},\qquad w=H_Bv,\qquad
\tau=-\frac{2B^TQw+\ell^Tw}{w^TQw},\qquad C=B+\tau w.
\]

完整狀態為 \(X_k=(P_{k-1},P_k)\)，則在每個定義良好的步驟上：

\[
\begin{bmatrix}P_{k+1}\\P_k\end{bmatrix}
=\underbrace{\begin{bmatrix}I+\tau_kH_k&-\tau_kH_k\\I&0\end{bmatrix}}_{M(X_k)}
\begin{bmatrix}P_k\\P_{k-1}\end{bmatrix}.
\]

注意 **\(M(X_k)\) 依賴撞擊狀態**，因此不能直接用 \(M(X_1)^n\) 代替整條撞擊鏈。程式不會為除法失效的步驟任意編造下一個點。

## 直接執行

需求：Python 3.10+；核心程式及測試不需第三方套件。程式可對大有限域進行**單條軌道反射**（超過 64 bits 的質數檢查僅為 probable-prime 測試）；`scan` 枚舉整個狀態圖，只適合小質數。

```bash
# 在專案根目錄執行
python -m unittest discover -s tests -v

# 搜尋 F_103 中三種曲線的週期圖
python -m conic_billiard scan --p 103 --kind all

# 只搜尋雙曲線
python -m conic_billiard scan --p 127 --kind hyperbola

# 雙曲線的雙點 2-cycle 示範
python -m conic_billiard trace --kind hyperbola --p 103 --A 5 0 --B 98 0 --steps 8

# 橢圓上經過 Jacobi 提升的玩具 DH（含 BSGS 破解）
python -m conic_billiard.toy_dh

# 反射逐步法、4x4 矩陣逐步法與 Jacobi 快速跳步計時
python -m conic_billiard.benchmark
```

`scan` 會列出 `max_period`、`largest_prime_factor_of_max_period` **以及** `largest_prime_factor_any_period`。最大循環不一定包含最大的質數階因子；後者對密碼學候選搜尋更有意義，但仍遠不足以評估安全性。

## 小域回歸結果

使用 \(a=5,b=3\)；結果由 `python -m conic_billiard scan` 算出，測試會重新核對。

| \(p\) | 曲線 | 仿射點數 | 最大週期 | 無法以仿射公式繼續的起始狀態 |
|---:|---|---:|---:|---:|
| 103 | 橢圓 | 104 | 60 | 0 |
| 103 | 雙曲線 | 102 | 60 | 196 |
| 103 | 拋物線 | 103 | 104 | 102 |
| 127 | 橢圓 | 128 | 126 | 0 |
| 127 | 雙曲線 | 126 | 74 | 244 |
| 127 | 拋物線 | 127 | 128 | 126 |

這些數字只屬於**指定的小質數、方程與仿射反射映射**。不能推論雙曲線普遍比橢圓好、或反之。對 \(p=103\)，兩條曲線的最大週期都是 60（最大質因子 5），但也都有週期 58 的軌道（最大質因子 29）。

## 專案組成

```
conic_billiard/
  core.py       # 有限域反射、原始 4x4 狀態矩陣、異常分類
  search.py     # 狀態圖、循環與質因子統計
  jacobi.py     # 橢圓型的已知虧格 1 群快速算法（對照基線）
  toy_dh.py     # 易破解的小域 ECDH 型示範
  benchmark.py  # 三種演算法的時間比較
  real_oval.py             # 實數封閉分支的物理第一撞擊點
  real_oval_experiment.py  # 反射選點、可逆性、數值穩定度
  real_oval_high_precision.py  # 可選 100 位數精度校驗
  __main__.py   # CLI

tests/          # 單元測試與已知數值回歸
docs/
  hyperbola.md      # 雙曲線與有限域定義域
  open_questions.md # 矩陣倍增、射影延拓、分類及密碼學檢驗
  security.md       # DH 作用、威脅模型及已知弱點
  references.md     # 文獻入口
```

## 兩件刻意沒有宣稱的事情

1. **沒有宣稱雙曲線已具有完整、安全的 ECDH 實作。** 雙曲線目前有原始反射與循環實驗；若要提升到虧格 1 群，還需針對具體不變層構造、驗證映射及群階。
2. **沒有宣稱已從原始 4×4 矩陣找到新的 \(O(\log n)\) 快速倍增。** `jacobi.py` 是已知群加法的比較基線，不是矩陣倍增的解答。這個問題在 `docs/open_questions.md` 被明確定義為研究任務。

## 研究資料

- [雙曲線完整數學定義與例外情況](docs/hyperbola.md)
- [矩陣倍增及其他待研究問題](docs/open_questions.md)
- [密碼學威脅模型與誤區](docs/security.md)
- [實數封閉分支：第一撞擊點與有限域的差異](docs/real_oval_first_hit.md)
- [參考資料](docs/references.md)

授權：MIT。可用於教學與探索，不可直接用來交換真正的金鑰。

## 實驗：有限域的「固定選正根」

對三次曲線反射後的二次方程，額外比較三種人為分支：固定採用正號（先以較小整數代表選平方根）、選取正的平衡代表值 `1 <= t <= (p-1)//2`、依撞擊點字典序選最小值。這些規則可以使部分狀態的下一點確定，但**不會自動保證可逆、幾何不變或有理交點存在**。

```bash
python -m conic_billiard.positive_branch_experiment
```

在 `p=101`、`E: y²=x³−x` 的 10,506 個不同的有序撞擊狀態中，固定使用「+」的規則只有 5,406 個定義良好的後繼、最長循環 10，且產生了 1,339 個入度大於 1 的狀態。固定選「+」的結果甚至可能因為將同一反射方向向量乘以非零純量而改變，所以它不是內在幾何的第一撞擊點定義。詳見 `results/positive_branch_p101.txt`。

## Finite-field group-order scan (ellipse and hyperbola)

We now screen group orders for both `epsilon=+1` (ellipse) and `epsilon=-1` (hyperbola). The lifted state space is a generalized Jacobi intersection whose associated elliptic curve has the same point count as `y^2 = x(x-epsilon)(x-m)` over the finite field. **The group order depends on `p`, the boundary, and the caustic/initial direction, not on `p` alone.** A promising DH toy subgroup has `#E(F_p)=4*r` with `r` a large prime; the actual reflected translation must also have a nontrivial `r`-component.

- [Exact group-order report](docs/group_order_search.md)
- [Order candidate CSV](results/conic_prime_orders.csv)
- [Order candidate JSON](results/conic_prime_orders.json)

Optional screening tool (requires `pip install sympy`):

```bash
python tools/conic_field_order_scan.py --primes 103 127 257 1009 1019 65537
```

For 255-/256-bit field primes, exact point counting requires optimized elliptic-curve routines (e.g. SageMath) and has **not** been completed by the included small/mid-size experiments:

```bash
sage -python tools/sage_large_prime_scan.py --prime curve25519 --trials 30
sage -python tools/sage_large_prime_scan.py --prime secp256k1 --trials 30
```

**Research only. No security claim; not suitable for production key exchange.**

## Large-prime sample verification

Four concrete 127/255/256-bit finite-field mirror-billiard samples, with exact related elliptic-curve orders and 64-step direct-reflection checks, are in [docs/large_prime_samples.md](docs/large_prime_samples.md) and [results/large_prime_samples.csv](results/large_prime_samples.csv). Two samples have certified or published prime-order subgroups (127-bit `20r`; 255-bit `8r`), while secp256k1- and P-256-*field* CM samples have composite `N/4` and remain unsuitable as selected DH parameters. Run `python tools/verify_large_prime_samples.py`.


## New: compact lifted-state encoding (32 bytes at 255 bits)

For the fixed-caustic generalized Jacobi model

\[c^2+\epsilon s^2=1,\quad d^2+m s^2=1,\]

one field parameter `t=s/(1+c)` and one root branch bit suffice to encode a full affine lifted state. Some projective exceptional points require a separate chart. A reserved `t=p` payload represents the two points with `c=-1` that lie outside the affine rational chart. All bytes are canonical, and invalid points/poles/square-root choices are rejected. This format works at 32 bytes for `p=2^255-19` and at a smaller size for toy fields; it is **not compatible with RFC 7748/X25519**.

```bash
python -m conic_billiard.compressed_dh_demo
python -m unittest discover -s tests -v
```

See [docs/compressed_points.md](docs/compressed_points.md) for the complete format and warnings. Encoding is not a substitute for subgroup validation or secure, constant-time arithmetic. This is an educational demo, **not deployable cryptography**.
