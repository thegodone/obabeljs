#!/bin/bash

rm obabel.js

cmake  \
-DCMAKE_TOOLCHAIN_FILE=//Users/mbp/Github/emsdk_portable1/emscripten/tag-1.34.6/cmake/Modules/Platform/Emscripten.cmake \
-DObabel_INCLUDE_DIR=/Users/mbp/Github/openbabel/include \
-DObabel_LIB_DIR=/Users/mbp/Github/build/src \
-DEMSCRIPTEN_BIN=/Users/mbp/Github/emsdk_portable1/emscripten/tag-1.34.6

make -j7


cat javascript/pre.js obabel.js javascript/post.js > dist/obabel.js

