# Phase 2: First-Principles Implementation (原理実装)
> **Action:** 作れなければ理解した扱いにするな．C言語の構造体，あるいは擬似コードで定義せよ．

## 1. The Struct (構造体定義)
- 最小の「状態」と「パラメータ」を `struct` に落とせ．(余計なメンバは入れるな)
```c
typedef struct {
    // 必要なメンバ変数は何か
    // double energy;
    // complex double state;
} science1_t;
```

## 2. The Algorithm (アルゴリズム)
- 入力を出力へ写す `Core Kernel` を書け．文法よりロジックを優先し，入出力の意味をコメントで固定せよ．
```c
void update(science1_t *sys) {
    // 入力: sys の現在状態
    // 出力: sys の更新後状態
    // ここに物理法則/ロジックを記述
}
```

## 3. Lecture Note (他者への説明)
- 実装の「肝」を30字以内で言い切れ．(条件分岐や例外ではなく，不変な中核だけ)
    -
