#! /bin/sh
# bootstrap file to be used when autogen.sh fails.
# consider using --foreign flat to automake
echo "Running aclocal..."
aclocal || exit 1
echo "Running autoheader..."
autoheader || exit 1
echo "Running autoconf..."
autoconf || exit 1
echo "Running automake..."
automake --foreign --add-missing --copy || exit 1
# automake-1.8 does not copy all the required files properly.
test -f config.sub || {
  echo "Your automake did not generate config.sub - this used to be an issue."
  echo "You may generate this file by running eg 'libtoolize -c'"
}
echo "Running configure $* ..."
./configure "$@"
