export PYTHONPATH=$HOME/codes/far/build-release:$HOME/codes/far/PHARE/pyphare:$PYTHONPATH
cd $HOME/codes/far/build-release
cmake -DCMAKE_CXX_FLAGS="-g3 -O3 -march=native -mtune=native -DPHARE_DIAG_DOUBLES=1" \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
      -DwithCaliper=OFF \
      -DtestMPI=OFF \
      -Dasan=OFF \
      -DdevMode=OFF \
      "$HOME"/codes/far/PHARE
make -j
