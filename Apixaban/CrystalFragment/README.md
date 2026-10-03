# Apixaban 結晶断片クラスター（案B）

## 目的
結晶構造から切り出した2分子（ダイマー）について、
1. 会合エネルギー ΔE_int = E(AB) − E(A) − E(B)
2. 1分子を剛体的に引き離したときのエネルギー変化（簡易な脱離曲線）
を、XTB2 / r2SCAN-3c の単一点計算で求める。

## 手順
1. Apixaban の結晶構造 CIF を入手し、このフォルダに置く
   - 対象: Form N-1（安定晶形）。CIF の入手元（CSD / 文献の補足情報など）は要確認
2. 必要なパッケージ: `pip install ase scipy numpy`
3. 入力ファイル生成（このフォルダで実行）
   ```
   python ..\..\scripts\make_crystal_fragments.py Apixaban_FormN1.cif --out . --top 3
   ```
   - `summary.csv`: 見つかった隣接ダイマーの一覧（重心間距離、最短距離、接触数）
   - `D1/ D2/ D3/`: 各ダイマーの入力ファイル
     - `*_complex_*`, `*_monoA_*`, `*_monoB_*`（XTB2 と r2SCAN-3c）→ ΔE_int 用
     - `*_pull_+X.XA_XTB2_SP.inp` → 引き離し曲線用
   - `run_all.bat`: 全ジョブを順に実行（orca が PATH にあること）
4. 実行後、ΔE_int = E(complex) − E(monoA) − E(monoB) を求める

## 注意
- 結晶構造のH原子位置は X 線由来だと不正確。絶対値を議論する前に、重原子固定でH原子だけ緩和することを検討
- 引き離しは剛体変位の単一点なので、あくまで概算（緩和・溶媒は含まない）
- 溶媒殻付きの評価は、ここで有望なダイマーを選んだあとの次段階（Explicit10の延長）
- ダイマー選びは接触数によるヒューリスティック。最終的な選択は summary.csv を見て判断すること
