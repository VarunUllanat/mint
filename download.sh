#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <host>"
  exit 1
fi

HOST="$1"
REMOTE_DIR="$(basename "$PWD")"

rsync -azvhP --existing "$HOST":"$REMOTE_DIR"/ .