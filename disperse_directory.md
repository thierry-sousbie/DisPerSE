**# 現在のディレクトリと研究方針**

このリポジトリ直下を、提示された `DisPerSE/` に相当するルートとして整理しています。

以下は現在の配置の概要です。後半の詳細ツリーも同じ実際の配置に揃えています。

```text

DisPerSE/（現在のリポジトリ直下）

├── 01_extraction/                 # DisPerSE strict/loose抽出

├── 02_label_generation/           # 共通関数・3クラスGT生成

├── 03_dataset/                    # 最終空間分割・パッチ保存・分割確認

├── 04_comparison/                 # 物理幅・FilFinder/DisPerSE・Herschel比較

├── 05_centerline_experiments/     # 共有ASPP予測の中心線可視化

├── archive_distance_and_width/    # 距離・幅・領域別の旧比較実装

└── archive_old_experiments/

├── 01_old_filament_comparison/

├── 02_old_dataset_split/

├── 03_old_prediction_visualization/

├── 04_old_evaluation/

├── 05_utilities/

└── 06_other_experiments/      # 点出力の確認・旧Orion Bスクリプト

```

**## 現方針との対応**

| 工程 | 現在の配置・対応コード | 比較・検証用の保存先 |

|---|---|---|

| FITS柱密度マップ → FilFinder/DisPerSE抽出 | `01_extraction/disperse_simajiri.py`、`disperse_simajiri_loose.py`。FilFinder抽出コードは別アップロード | `archive_distance_and_width/` に距離・幅・領域別実装 |

| 3クラスGT生成（0:背景、1:Class1、2:Class2、255:無効） | `02_label_generation/generate_and_visualize_3class_map.py` と同階層の `label3class_common.py` | `04_comparison/`、`archive_old_experiments/01_old_filament_comparison/` で抽出手法を比較 |

| 空間分割 | `03_dataset/make_regions_connected_no_leakage_v5.py` | `archive_old_experiments/02_old_dataset_split/` の初期分割・斜め帯・連結成分別実装 |

| 入力・GTのパッチ化（split_v2） | `03_dataset/save_final_patch_mask_pairs_chunked.py`。分割確認は同階層の `visualize_saved_split_*.py` | 旧一体型パッチ生成は `archive_old_experiments/02_old_dataset_split/` |

| 初期：Class1とClass2の独立2値分類 | 学習コードは別ZIP | このリポジトリのみではモデル比較を実行できない |

| 共有特徴抽出器＋2-headモデル | 学習コードは別ZIP | 旧Validation PRは `archive_old_experiments/04_old_evaluation/plot_pr_shared_aspp_class1_class2_validation.py` |

| checkpoint整合・損失・評価処理の修正 | 修正版は別ZIPで確認 | 使用config・checkpoint・コード版の対応を記録 |

| HRF事前学習＋共有ASPP／one-hot CE＋Dice | 学習・損失実装は別ZIP | `05_centerline_experiments/visualize_shared_aspp_class1_class2_connected_lines.py` は中心線表示の検証用 |

| PR曲線・F1評価・元FITS座標で可視化 | 現方針の評価コードは別ZIP | 旧PRは `04_old_evaluation/`、旧座標合わせ・Class1/Class2表示は `03_old_prediction_visualization/`（いずれも `archive_old_experiments/` 配下） |

| 後続：Head1=Class1、Head2=Class1∪Class2（Union）の比較 | 2026年10月下旬の別ZIP。現フォルダには未収録 | 同一split・GT・評価条件で従来2-headと比較する |

**# DisPerSE.zip：ディレクトリ構造と研究処理の整理**

**# 現在のディレクトリと研究方針**

```text

DisPerSE/

├── 01_extraction/

│   ├── disperse_simajiri.py  # DisPerSE strict：FITS→mse→VTK→線・マスク・画像

│   └── disperse_simajiri_loose.py  # DisPerSE loose：検出しきい値を緩めた版

├── 02_label_generation/

│   ├── generate_and_visualize_3class_map.py  # strict/looseを統合して3クラスFITS・PNGを生成

│   └── label3class_common.py  # 距離・画素スケール・膨張・短い線の除去など共通関数

├── 03_dataset/

│   ├── make_regions_connected_no_leakage_v5.py  # 分子雲ごとの空間分割を確定、漏洩防止

│   ├── save_final_patch_mask_pairs_chunked.py  # 保存済みmanifestから入力/GTパッチをPNG化

│   ├── visualize_saved_split_boundaries_original_colors.py  # FITS上に分割境界を描く

│   └── visualize_saved_split_islands_png.py  # 分割パッチの飛び地・孤立を可視化

├── 04_comparison/

│   ├── compare_disperse_filfinder_fixed_red.py  # 両手法の一致・不一致を色分けした修正版

│   ├── compare_skeleton.py  # Herschel公開スケルトンとの視覚的比較

│   ├── compare_vs_herschel.py  # Herschel公開結果との比較

│   └── disperse_filfinder_comparison_physical.py  # 物理幅マスクを使う比較

├── 05_centerline_experiments/

│   └── visualize_shared_aspp_class1_class2_connected_lines.py  # 予測領域を線状にしてGTと重ねる

├── archive_distance_and_width/

│   ├── disperse_aquila.py  # AquilaのDisPerSE抽出を試す

│   ├── disperse_aquila_physical.py  # 物理幅を考慮したDisPerSE抽出の試作

│   ├── disperse_filfinder_comparison.py  # 両手法の中心線・マスク比較

│   ├── disperse_filfinder_simajiri.py  # 2手法の一致・不一致を色分け

│   └── disperse_simajiri_orionB_fixed.py  # Orion Bの大きい画像を縮小処理してDisPerSEを実行

├── archive_old_experiments/

│   ├── 01_old_filament_comparison/

│   │   └── disPerse_filfiinder_comparison.py  # FilFinder対DisPerSEの初期比較

│   ├── 02_old_dataset_split/

│   │   ├── make_region_3class_patches_chunked.py  # 分割・正規化・256pxパッチ保存を一体化

│   │   ├── make_region_3class_patches_chunked_no_hdu.py  # 連結成分を考慮した分割の別実装

│   │   ├── split_aquila_diagonal.py  # 斜め帯での分割とフィラメント越境対策

│   │   ├── train_val_test.py  # Aquilaの空間・クラス比を考慮した初期分割

│   │   └── train_val_test_png.py  # 分割結果の可視化のみ

│   ├── 03_old_prediction_visualization/

│   │   ├── visualize_candidate_test_contour.py  # Class2のTest輪郭表示

│   │   ├── visualize_candidate_test_from_patch_masks.py  # Class2のTest GT/予測を再配置

│   │   ├── visualize_confirmed_correct_coordinates.py  # 元FITS座標系に合わせてClass1予測を表示

│   │   ├── visualize_confirmed_test_contour.py  # Class1のTest輪郭表示

│   │   └── visualize_confirmed_test_from_patch_masks.py  # Class1のTest GT/予測を再配置

│   ├── 04_old_evaluation/

│   │   └── plot_pr_shared_aspp_class1_class2_validation.py  # 共有ASPPモデルのval PR曲線

│   ├── 05_utilities/

│   │   ├── npy_open.py  # NumPy配列の中身を確認

│   │   └── seg.py  # 3クラスセグメンテーションの説明用図

│   └── 06_other_experiments/

│       ├── disperse_point.py  # DisPerSEの点・出力の可視化確認

│       └── disperse_simajiri_orionB.py  # 空ファイル（0 byte）

├── .gitattributes

├── .gitignore

├── cleanup

├── CMakeLists.txt

├── disperse_directory.md

├── readme (2).md

└── readme.md

```

**## 2. FITSから学習パッチができるまで**

```text

元の柱密度FITS

├─ FilFinder / DisPerSE のstrict・loose結果

│       └─ generate_and_visualize_3class_map.py

│            └─ 3クラス暫定ラベルFITS

│

└─ 元の柱密度FITS（学習画像）

│

▼

make_regions_connected_no_leakage_v5.py

└─ train / val / test の空間分割を決定

（フィラメント連結成分の分割またぎを考慮）

│

▼

metadata/<region>/patch_manifest.csv

└─ 各パッチのsplitと切り出し座標を保存

│

▼

save_final_patch_mask_pairs_chunked.py

├─ FITSの柱密度をlog10変換・正規化して8bit画像化

├─ 同じ座標から暫定ラベルを切り出す

└─ 256×256pixelの入力PNG・GT PNGを保存

│

▼

split_v2/{train,val,test}/{in,out}/<region>/...

│

▼

MMSegmentationで学習・評価

``

**### パッチの保存構造（スクリプトが想定する出力例）**

```text

split_v2/

├── metadata/

│   └── aquila/

│       └── patch_manifest.csv                # split・座標の管理表

├── train/

│   ├── in/aquila/chunk_0000/y01000_x02000_in.png  # 正規化柱密度画像

│   └── out/aquila/chunk_0000/y01000_x02000_out.png # 3クラスGT

├── val/

│   ├── in/...

│   └── out/...

└── test/

├── in/...

└── out/...

```

**## 3. スクリプト保存日時による研究の流れ**

| ZIP内更新日（概略） | 代表スクリプト | 読み取れる研究内容 |

|---|---|---|

| 2026-07-03〜08 | `disperse_aquila*`, `disPerse_filfiinder_comparison.py`, `disperse_filfinder_comparison_physical.py` | AquilaのDisPerSE抽出、FilFinderとの中心線・物理幅マスク比較 |

| 2026-07-17〜24 | `disperse_simajiri*`, `compare_*` | 他領域対応、境界ノイズ処理、公開スケルトンとの比較・可視化 |

| 2026-07-26〜29 | `train_val_test*`, `label3class_common.py`, `generate_and_visualize_3class_map.py`, `split_aquila_diagonal.py` | strict/looseによる3クラスラベル生成、初期の空間分割とフィラメント越境対策 |

| 2026-08-06 | `make_region_3class_patches_chunked*.py` | 分子雲ごとの分割、log10正規化、256px学習パッチ化の検討 |

| 2026-08-27〜31 | `visualize_confirmed_*`, `visualize_candidate_*` | Testパッチの元画像座標への再配置、GTと予測のずれの確認 |

| 2026-09-01 | `make_regions_connected_no_leakage_v5.py`, `save_final_patch_mask_pairs_chunked.py`, `visualize_saved_split_*` | 連結フィラメントを意識した最終分割と、manifestに基づくパッチ保存・分割検査 |

| 2026-10-06 | `plot_pr_shared_aspp_class1_class2_validation.py`, `visualize_shared_aspp_class1_class2_connected_lines.py` | 共有ASPP型2-headのValidation PR曲線と予測中心線風の可視化 |

**## 4. 再現実験前の重要な確認点**

1. ****距離・幅の整合性****：分子雲ごとの距離、`WIDTH_PC`、膨張半径・カーネル、DisPerSEの `thickness=3` が教師ラベル生成時にどう作用するか確認。

2. ****元の3クラスラベルと最終パッチGTの一致****：Class1＝strict和集合、Class2＝loose和集合−Class1、255＝ignoreの扱いを確認。

3. ****分割の再現性****：`patch_manifest.csv` が実際に使ったものか、連結成分が複数splitに跨がっていないか、train/val/testの重複がないか確認。

4. ****正規化****：Trainだけから正規化範囲を推定しているか、validation/testが同じ範囲を利用しているか確認。

5. ****細線化可視化と領域評価の区別****：線として表示した画像の見た目と、領域ベースのF1/PR曲線を混同しない。

6. ****使用したコード版の確定****：同じ役割の複数スクリプトがあるため、実際に使用した版とチェックポイントを記録する。