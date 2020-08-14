#!/bin/bash

find src -name "*.cpp" -or -name "*.h" -exec clang-format -i {} \;
