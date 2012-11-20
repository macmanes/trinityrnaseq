if [ -z "$nCPUs" ]; then
    nCPUs=$(grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu)
fi
pref=$(basename $0 .sh)
DIR=../bin
JF=$DIR/jellyfish

check () {
    cut -d\  -f 2 $1 | xargs md5sum | sort | diff -w $1 -
}

if [ -n "$DEBUG" ]; then
    set -x;
fi

set -e
