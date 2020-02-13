using Documenter
makedocs(;
    sitename="Juliaで学ぶ量子力学",
    pages = [
        "Index"=>"index.md",
        "時間に依存しないシュレーディンガー方程式の解"=>[
            "ポテンシャルがない場合1次元シュレーディンガー方程式を解き、その後数値的に解いてみる" => "chapter1/01.md",
            "ポテンシャルがある場合の1Dシュレーディンガー方程式の解" => "chapter1/02.md",
            #"About Goma-chan" => "chapter1/goma.md",
        ],
#        "Chapter2"=>"chapter2/azarashi.md",
    ]
)