#!/bin/bash

#switch to the folder of this script
DESTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DESTDIR
BINDIR=$DESTDIR/bin
DOCDIR=$DESTDIR/doc
MANDIR=$DESTDIR/man
rm -rf $BINDIR
rm -rf $DOCDIR
rm -rf $MANDIR

FILE=`ls -d collectl*| tail -1`
tar xzf ${FILE}
INSTALLDIR=`find . -maxdepth 1 -type d -iname "*collect*" | tail -1`

echo "DESTDIR=$DESTDIR"
echo "FILE=$FILE"
echo "INSTALLDIR=$INSTALLDIR"

mkdir -p $BINDIR
mkdir -p $DOCDIR
mkdir -p $MANDIR

cp make_data_files.sh    $BINDIR
cp plot.sh               $BINDIR
cp timetable.sh          $BINDIR
cd ${INSTALLDIR}
cp collectl.pl           $BINDIR/collectl
cp collectl.conf         $BINDIR
cp -r man1               $MANDIR

cp docs/*                $DOCDIR
cp GPL ARTISTIC COPYING  $DOCDIR
cp RELEASE-collectl      $DOCDIR

cp UNINSTALL             $BINDIR
cp *.ph                  $BINDIR
cp envrules.std          $BINDIR
cp client.pl             $BINDIR
cp col2tlviz.pl          $BINDIR

chmod 444 $BINDIR/collectl.conf
chmod 755 $BINDIR/collectl
chmod 444 $DOCDIR/ARTISTIC $DOCDIR/COPYING $DOCDIR/GPL
chmod 755 $BINDIR/*.pl
chmod 444 $BINDIR/*.ph

cd ${DESTDIR}
rm -rf ${INSTALLDIR}


