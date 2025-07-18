c++ バージョンと Rust バージョンがある．c++ バージョンは make で，Rust バージョンは cargo run --release でコンパイル可能である．

Rust install: https://www.rust-lang.org/ja/tools/install

frinZ.py の拡張・補助ツール．frinZ.py は精密にフリンジの位置を推定できないので，それを補うのが frinZsearch である．     
frinZsearch が推定した delay と rate を用いて frinZ.py を実行することで，遅延較正・位相較正が可能となる．   


コンパイル（make --> make install）は ubuntu 24.04 (GCC 13)と Rocky linux 8.10 (GCC 8) で試しました．

FFTW を用いているので，-lfftw3f が 適切な箇所に保存もしくはリンクがないと make できません．-lfftw3f で make できないときは， g++ -print-search-dirs で g++ が利用するライブラリーの PATH を確認してください．適当に /usr/lib にリンクを貼れば大丈夫だと思います．

\#include \<filesystem\> を用いているので，g++ のバージョンによっては -lstdc++fs が必要（GCC 8.x 以前？）なので，一応 LDFLAGS には追加しています．

