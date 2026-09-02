# -*- coding: utf-8 -*-
"""
solvent_swap.py

「既存の明示的溶媒和クラスター（例: Apixaban + 4 H2O）」の座標を土台にして、
水素結合の方向（O-H --- 溶質の極性部位）を保ったまま、別の溶媒分子
（例: メタノール, O-CH3 / DMSO, S(=O)(CH3)2）に置き換えた初期構造を生成するスクリプト。

背景・設計方針
--------------
- 水分子は「極性部位近くでクラッシュを避けて配置済み」という前提のもと、
  各水分子の O 原子位置と、溶質側を向いている O-H 結合方向（水素結合ドナー方向）
  はそのまま維持する。
- メタノール: 溶質と反対を向いているもう一方の H だけを、CH3基に置き換える。
- DMSO: 水のOをDMSOのS=Oの酸素として流用（配位方向は水と同じ役割）し、
  溶質と反対方向にS原子を配置、Sにぶら下がる2つのメチル基は実測に近い
  O-S-C角(~106.7°)・C-S-C角(~97.4°)のピラミッド型構造で追加する。
- どちらも、追加した原子が既存原子（溶質 or 他の溶媒分子）とクラッシュ
  （原子間距離が近すぎる）していないかを自動チェックし、必要なら結合長を
  少しずつ伸ばして回避する（ORCAのOptimizerが後で微調整する前提の
  「まず壊れない初期構造」を作ることが目的）。

使い方（例）
------------
    python solvent_swap.py \
        --input  Apixaban_Explicit10H2O_repositioned_initial.xyz \
        --output Apixaban_Explicit10DMSO_repositioned.xyz \
        --n-solvent 10 \
        --template dmso \
        --bond 1.53

--n-solvent には、xyzファイルの「末尾から」何分子が明示的溶媒(水)かを指定する。
（このプロジェクトの命名規則上、Explicit4H2O等は常に O,H,H の3行×N分子が
  ファイル末尾に並んでいる）

このスクリプトは標準的なxyz座標の幾何操作のみを行う。ORCAの.inpファイル生成は
含まない（.inpは計算条件〈汎関数・基底関数・SMD溶媒名など〉が都度変わるため、
座標だけをここで確定させ、.inpのヘッダーはClaude側で軽く仕上げる想定）。

更新履歴
--------
- v2: クラッシュ判定のバグを修正。以前のバージョンは、追加した原子（C, S）を
  「自分自身が結合しているO原子」に対してもクラッシュ判定してしまい、
  共有結合の距離（1.4〜1.5Å程度）をvdW半径の和（2.4Å前後）と比較して
  常に「クラッシュ」と誤判定していた（毎回結合長が伸びてしまう不具合）。
  現在は、新しく置換基を作る際、その"親"であるOおよび保持したHをクラッシュ
  判定の対象から除外し、他の原子（溶質・他の溶媒分子）とのみ比較するように修正済み。
"""

import argparse
import numpy as np

# 元素ごとのvan der Waals半径 (Å) — クラッシュ判定に使用（簡易版）
VDW_RADII = {
    "H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
}

# 溶媒テンプレート定義:
#   置換するH原子の代わりに追加する置換基の情報。
#   "anchor_bond": 溶質と反対側のH原子だったO位置から新原子までの既定結合長 (Å)
#   "builder": 骨格原子(O位置, 保持するH-結合方向ベクトル, 除去したH方向ベクトル)
#              から置換基の原子群を作る関数
SOLVENT_TEMPLATES = {}


def register_template(name):
    def deco(fn):
        SOLVENT_TEMPLATES[name] = fn
        return fn
    return deco


def _perp_basis(axis):
    """axisに垂直な正規直交基底2本を返す"""
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, axis)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    perp1 = np.cross(axis, tmp)
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(axis, perp1)
    return perp1, perp2


def _attach_ch3(anchor_pos, bond_dir, ch_bond=1.09, phi0=0.0):
    """
    anchor_pos（C原子の位置）から、bond_dir方向を「元になった結合の延長方向」として
    3つのHをtetrahedral配置で追加する。methanolのCH3構築ロジックを汎用化したもの。
    phi0: 3本のHの方位角オフセット（他の置換基との干渉を避けたい場合に使う）
    戻り値: [("H", xyz), ("H", xyz), ("H", xyz)]
    """
    perp1, perp2 = _perp_basis(bond_dir)
    tetra_angle = np.radians(109.5)
    atoms = []
    for k in range(3):
        phi = np.radians(120.0 * k) + phi0
        direction = (
            np.cos(tetra_angle) * bond_dir
            + np.sin(tetra_angle) * (np.cos(phi) * perp1 + np.sin(phi) * perp2)
        )
        direction /= np.linalg.norm(direction)
        atoms.append(("H", anchor_pos + direction * ch_bond))
    return atoms


@register_template("methanol")
def build_methanol(o_pos, away_dir, anchor_bond=1.43):
    """
    O位置(o_pos)から、溶質と反対方向(away_dir, 単位ベクトル)にO-C結合を伸ばし、
    メタノールのCH3基を構築する。
    O-C結合長の既定値は1.43 Å（メタノールの実測値近傍）。
    Cにぶら下がる3つのHは、O-C軸に対して交互(スタガード)になるよう
    tetrahedral配置で追加する。
    戻り値: [(element, xyz), ...] の追加原子リスト（Oは含まない、Cと3Hのみ）
    """
    c_pos = o_pos + away_dir * anchor_bond
    atoms = [("C", c_pos)]
    atoms.extend(_attach_ch3(c_pos, away_dir))
    return atoms


@register_template("dmso")
def build_dmso(o_pos, away_dir, anchor_bond=1.53):
    """
    O位置(o_pos)を保持し、水のH-結合方向と同じ役割（溶質への配位方向）を
    DMSOのS=O酸素にも担わせる。S原子はaway_dir方向に伸ばし(既定1.53 Å = S=O結合長)、
    Sにぶら下がる2つのメチル基は、実測に近いO-S-C角(~106.7°)・C-S-C角(~97.4°)を
    再現するように配置する（Sはピラミッド型）。
    anchor_bond引数はここではS=O結合長として使う（クラッシュ回避時に伸長される）。
    戻り値: [(element, xyz), ...] の追加原子リスト（Oは含まない、S,C,C,6Hの9原子）
    """
    s_pos = o_pos + away_dir * anchor_bond

    perp1, perp2 = _perp_basis(away_dir)

    # S→O方向 (-away_dir) を基準に、O-S-C角 ≈106.7°だけ開いた2方向にCを置く。
    # 2本のC-S結合のなす角(C-S-C)が実測値97.4°になるよう、方位角差を球面幾何から逆算。
    o_s_c_angle = np.radians(106.7)
    cos_csc = np.cos(np.radians(97.4))
    cos_osc = np.cos(o_s_c_angle)
    sin_osc = np.sin(o_s_c_angle)
    cos_delta = (cos_csc - cos_osc ** 2) / (sin_osc ** 2)
    cos_delta = max(-1.0, min(1.0, cos_delta))
    delta_phi = np.arccos(cos_delta)

    c_bond = 1.80  # S-C結合長 (Å)
    atoms = [("S", s_pos)]
    c_positions = []
    for sign in (+1, -1):
        phi = sign * delta_phi / 2.0
        direction = (
            np.cos(o_s_c_angle) * (-away_dir)
            + np.sin(o_s_c_angle) * (np.cos(phi) * perp1 + np.sin(phi) * perp2)
        )
        direction /= np.linalg.norm(direction)
        c_pos = s_pos + direction * c_bond
        c_positions.append(c_pos)
        atoms.append(("C", c_pos))

    # 各メチル基にH3個ずつ（C-S結合方向を軸にtetrahedral配置）
    for idx, c_pos in enumerate(c_positions):
        c_dir = c_pos - s_pos
        c_dir /= np.linalg.norm(c_dir)
        atoms.extend(_attach_ch3(c_pos, c_dir, phi0=np.radians(60.0 * idx)))

    return atoms


def read_xyz(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    n = int(lines[0].strip())
    comment = lines[1]
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        el = parts[0]
        xyz = np.array([float(v) for v in parts[1:4]])
        atoms.append((el, xyz))
    return comment, atoms


def write_xyz(path, comment, atoms):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for el, xyz in atoms:
            f.write(f"{el:<2s}  {xyz[0]:14.6f}  {xyz[1]:14.6f}  {xyz[2]:14.6f}\n")


def has_clash(new_atoms, existing_atoms, factor=0.75):
    """
    new_atoms: このステップで追加した原子 [(el, xyz), ...]
    existing_atoms: 比較対象の原子群（自分自身が直接結合しているO・保持したHは
                     呼び出し側で除外済みであること — 直接結合はvdW半径の和より
                     短くて当然なので、含めると常にクラッシュ判定されてしまう）
    factor: vdW半径の和に掛ける係数。小さいほど判定が緩くなる。
            0.75前後は「かなり接近しているが計算が破綻するほどではない」
            ラインの目安（会話ログでの2.58Å伸長の例に近い挙動）。
    """
    for el_new, xyz_new in new_atoms:
        r_new = VDW_RADII.get(el_new, 1.5)
        for el_ex, xyz_ex in existing_atoms:
            r_ex = VDW_RADII.get(el_ex, 1.5)
            dist = np.linalg.norm(xyz_new - xyz_ex)
            if dist < (r_new + r_ex) * factor:
                return True
    return False


def swap_solvent(input_path, output_path, n_solvent, template_name,
                  base_bond=1.43, bond_step=0.05, max_bond=3.0):
    comment, atoms = read_xyz(input_path)

    n_total = len(atoms)
    n_water_atoms = n_solvent * 3  # O,H,H ×n_solvent
    solute_atoms = atoms[: n_total - n_water_atoms]
    water_atoms = atoms[n_total - n_water_atoms:]

    builder = SOLVENT_TEMPLATES[template_name]

    new_atoms = list(solute_atoms)  # ここに確定した原子を積み上げていく
    report = []

    for i in range(n_solvent):
        o_el, o_pos = water_atoms[i * 3 + 0]
        h1_el, h1_pos = water_atoms[i * 3 + 1]
        h2_el, h2_pos = water_atoms[i * 3 + 2]
        assert o_el == "O"

        # 溶質(solute_atoms)の中で最も近い重原子(H以外)を探し、
        # そちらを向いているHを「保持するH（水素結合ドナー方向）」とみなす
        heavy_solute = [(el, xyz) for el, xyz in solute_atoms if el != "H"]

        def min_dist_to_solute(p):
            return min(np.linalg.norm(p - hp) for _, hp in heavy_solute)

        # H自身ではなく「Hの外側延長線上」で判定すると安定するため、
        # 単純にH原子自体の溶質重原子への最短距離で判定する
        if min_dist_to_solute(h1_pos) <= min_dist_to_solute(h2_pos):
            keep_h_el, keep_h_pos = h1_el, h1_pos
            away_h_el, away_h_pos = h2_el, h2_pos
        else:
            keep_h_el, keep_h_pos = h2_el, h2_pos
            away_h_el, away_h_pos = h1_el, h1_pos

        away_dir = away_h_pos - o_pos
        away_dir /= np.linalg.norm(away_dir)

        # クラッシュが解消するまで結合長を少しずつ伸ばす
        bond = base_bond
        while True:
            substituent = builder(o_pos, away_dir, anchor_bond=bond)
            # 自分自身のO/保持したH（=直接結合している相手）はクラッシュ判定から除外する。
            # 除外しないと「新しく作った共有結合」自体を誤ってクラッシュと判定してしまう。
            existing_for_check = new_atoms
            if not has_clash(substituent, existing_for_check):
                break
            bond += bond_step
            if bond > max_bond:
                report.append(
                    f"  [警告] 溶媒分子 #{i+1}: {max_bond}Åまで伸ばしてもクラッシュ回避不可。"
                    f" 手動確認が必要です。"
                )
                break

        if bond > base_bond:
            report.append(
                f"  溶媒分子 #{i+1}: 既定結合長 {base_bond:.2f}Å ではクラッシュしたため "
                f"{bond:.2f}Å に伸長（ORCA Optimizerが後で修正する想定）"
            )

        new_atoms.append((o_el, o_pos))
        new_atoms.append((keep_h_el, keep_h_pos))
        new_atoms.extend(substituent)

    write_xyz(output_path, comment + f" | solvent swapped to {template_name}", new_atoms)

    print(f"完了: {output_path} を書き出しました（原子数: {len(new_atoms)}）")
    if report:
        print("--- 注意事項 ---")
        for line in report:
            print(line)
    else:
        print("クラッシュなし。すべて既定結合長で配置できました。")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="土台となる既存xyzファイル（末尾がO,H,H×N分子の並び）")
    parser.add_argument("--output", required=True, help="出力先xyzファイル")
    parser.add_argument("--n-solvent", type=int, required=True, help="末尾から何分子が明示的溶媒(水)か")
    parser.add_argument("--template", default="methanol", choices=list(SOLVENT_TEMPLATES.keys()))
    parser.add_argument("--bond", type=float, default=1.43, help="既定の O-X 結合長 (Å)（DMSOはS=O結合長として使われる。目安1.53）")
    args = parser.parse_args()

    swap_solvent(args.input, args.output, args.n_solvent, args.template, base_bond=args.bond)


if __name__ == "__main__":
    main()
