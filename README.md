コンパイルは ubuntu 24.04 と Rocky linux 8.10 で試しました．

FFTW を用いているので，-lfftw3f が 適切な箇所に保存もしくはリンクがないと make できない．-lfftw3f で make できないときは， g++ -print-search-dirs で g++ が利用するライブラリーの PATH を確認してください．適当に /usr/lib にでもリンクを貼れば大丈夫だと思います．

\#include <filesystem> を用いているので，g++ のバージョンによっては -lstdc++fs が必要（GCC 8.x 以前？）なので，一応 LDFLAGS には追加している．

