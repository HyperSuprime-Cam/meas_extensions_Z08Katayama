# This config file can be used to enable the plugin algorithm when running
# the processCcd.py, scProcessCcd.py, and processCoadd.py scripts; use
# the "-C" option followed by the path to this file.
import lsst.meas.extensions.shapeZ08Katayama
root.measurement.algorithms.names |= ["shape.z08"]
