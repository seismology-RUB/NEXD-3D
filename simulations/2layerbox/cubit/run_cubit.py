#!python
#!/usr/bin/env python

import cubit
import cubit2dg3d
import os
import sys
from math import *
from numpy import *


reload(cubit2dg3d)

cubit.init([""])
command = 'open "2layerbox_PML.cub"'
cubit.cmd(command)

cubit2dg3d.mesh()
