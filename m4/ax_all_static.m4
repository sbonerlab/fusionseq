AC_DEFUN([AX_ALL_STATIC], [
	AC_ARG_ENABLE(all-static,
		       AC_HELP_STRING([--enable-all-static], 
				      [Build statically-linked-binaries]),
		       CC_ALL_STATIC=yes,
		       CC_ALL_STATIC=no)

	if test "$CC_ALL_STATIC" = yes; then
		OTHERLDFLAGS="-all-static"
		ROOTLDFLAGS="-L$ROOTLIBDIR -lcurses -lRoot -lpcre"
	else
		OTHERLDFLAGS=""
		ROOTLDFLAGS="-L$ROOTLIBDIR $ROOTGLIBS $ROOTLIBS -ldl"
	fi
	
	AC_SUBST([OTHERLDFLAGS])
	AC_SUBST([ROOTLDFLAGS])
])
