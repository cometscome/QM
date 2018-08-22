# Juliaで学ぶ量子力学  
In Japanese. Juliaで学ぶ量子力学  
最初は普通にシュレーディンガー方程式を解いて、そのあと、演算子表示を使って色々なやり方で解く予定。数値的に解く時には、Juliaを使用。  
随時更新予定。更新ペースはのんびり。

ファイルが重くなってしまったので、プレビュー用にPDFファイルを用意してみた。

こちらから見ると見やすいかもしれない。
http://nbviewer.jupyter.org/github/cometscome/QM/tree/master/

## 目次
### 時間に依存しないシュレーディンガー方程式の解
- 01 ポテンシャルがない場合1次元シュレーディンガー方程式を解き、その後数値的に解いてみる。
- 02 ポテンシャルがある場合の1次元シュレーディンガー方程式を数値的に解いてみる。
- 03 波数表示で解いてみる。ガウス関数形ポテンシャルのある問題
- 04 二次元シュレーディンガー方程式の解、平面波基底の解とベッセル関数基底の解。ポテンシャルがない場合のエネルギー固有値の比較
- 05 二次元シュレーディンガー方程式の解、平面波基底の解とベッセル関数基底の解。波動関数の比較。
- 06 二次元シュレーディンガー方程式を差分化して解いた時にエルミート行列になっていない問題について考察

#### Julia言語とは
Juliaは記述が簡単で高速な言語。行列の対角化から特殊関数まで、物理で用いる様々な計算を手軽に実行することができる。  
試す場合には  
https://www.juliabox.com/  
を使うとブラウザ上でこれらのノートブックを試すことができる。ログインして、ファイルをアップロードすれば、実行してみることが（もちろんいじることも）できる。

# Julia 1.0.0に順次対応中
