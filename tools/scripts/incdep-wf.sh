#!/bin/sh
# this is <incdep>
# ----------------------------------------------------------------------------
# $Id: incdep,v 1.2 2002/09/20 15:00:33 forbrig Exp $
# 
# Copyright (c) 2002 by Thomas Forbriger (IMG Frankfurt) 
# 
# look for fortran include dependencies
# 
# REVISIONS and CHANGES 
#    20/09/2002   V1.0   Thomas Forbriger
# 
# ============================================================================
#
for dir in . $@
do
for d in $dir/*.f
do
  fina=`basename $d .f`.o
  depna=`cat $d | grep '^	include' | cut -f 2 -d \' | sort | uniq`
  echo $fina: $depna
done
for d in $dir/*.f.m4
do
  fina=`basename $d .f.m4`.o
  depna=`cat $d | grep '^	include' | cut -f 2 -d \' | sort | uniq`
  echo $fina: $depna
done
for d in $dir/*.f90
do
  fina=`basename $d .f90`.o
  depna=`cat $d | grep '^	use' | cut -f 2 -d ' ' | cut -f 1 -d ',' | sort | uniq`
  for name in $depna
  do
	  echo $fina: $name.o
  done
done
done
# ----- END OF incdep ----- 
