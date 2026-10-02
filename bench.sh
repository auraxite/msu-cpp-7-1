#!/bin/bash
# Usage: ./bench.sh [file.tsv ...]  (default: small and large examples)

cd "$(dirname "$0")"
g++ -O2 -std=c++17 -o run main.cpp || exit 1
py=python3; "$py" -c "" 2>/dev/null || py=python
[ $# -eq 0 ] && set -- small_example.tsv large_example.tsv
for f in "$@"; do
  b="${f%.tsv}.bin"
  o="${f%.tsv}.out.tsv"
  t0=$(date +%s%N); ./run -s -i "$f" -o "$b"
  t1=$(date +%s%N); ./run -d -i "$b" -o "$o"
  t2=$(date +%s%N)
  ok=$("$py" -c "import sys
e = lambda p: sorted((min(a, b), max(a, b), w) for a, b, w in
                     (l.split('	') for l in open(p).read().splitlines()))
print('OK' if e(sys.argv[1]) == e(sys.argv[2]) else 'FAIL')" "$f" "$o")
  printf '%-20s %10d bytes  ser %5d ms  deser %5d ms  %s\n' \
    "$f" "$(wc -c < "$b")" \
    $(( (t1 - t0) / 1000000 )) $(( (t2 - t1) / 1000000 )) "$ok"
  rm -f "$o"
done
