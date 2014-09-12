#! /bin/sh

case `uname` in 
Darwin*) 
  alias cmd_libtoolize=glibtoolize ;;
*) 
  alias cmd_libtoolize=libtoolize ;;
esac

bootstrap() {
  cmd_libtoolize -c && \
  aclocal -I config && \
  automake -a -c && \
  autoconf
}

bootstrap
