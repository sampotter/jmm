#!/bin/bash

rm gtags.files

find src -name *.c >> gtags.files
find src -name *.cpp >> gtags.files
find src -name *.h >> gtags.files
find examples -name *.c >> gtags.files
find examples -name *.cpp >> gtags.files
find examples -name *.h >> gtags.files
find test -name *.c >> gtags.files
find test -name *.cpp >> gtags.files
find test -name *.h >> gtags.files
