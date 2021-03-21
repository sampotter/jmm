#!/usr/bin/env sh

OFF_PATH=$1

echo $OFF_PATH

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sed -i '/^#/d' $OFF_PATH
elif [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' '/^#/d' $OFF_PATH
else
    echo "Unknown platform: $OSTYPE" 1>&2
    exit 1
fi