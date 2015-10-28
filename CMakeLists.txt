set(obabeljs_SRC ${CMAKE_CURRENT_SOURCE_DIR}/obabeljs.cpp) 
set(obabeljs_H ${CMAKE_CURRENT_SOURCE_DIR}/obabeljs.h)


add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/obabel.js
        COMMAND ${EMSCRIPTEN_BIN}/em++ --bind -o obabel.js -s MODULARIZE=1 -O2 -s EXPORT_NAME='"'obabel'"'  ${obabeljs_SRC}  ${CMAKE_CXX_FLAGS}  -I${OBABEL_INCLUDE_DIR} 
        ${OBABEL_LIB_DIR}/libopenbabel.a --memory-init-file 0 -s DISABLE_EXCEPTION_CATCHING=0 -s TOTAL_MEMORY=16777216 -s ALLOW_MEMORY_GROWTH=1)
add_custom_target(obabel ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/obabel.js)

