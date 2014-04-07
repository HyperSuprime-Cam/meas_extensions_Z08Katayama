// -*- lsst-c++ -*-
%define shapeZ08Katayama_DOCSTRING
"
Access to the C++ classes from the meas_extensions_shapeZ08Katayama library
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.extensions.shapeZ08Katayama",
        docstring=shapeZ08Katayama_DOCSTRING) shapeZ08KatayamaLib

%{
#include "lsst/pex/logging.h"
#include "lsst/meas/algorithms.h"
#include "lsst/meas/extensions/shapeZ08Katayama.h"
%}

%include "lsst/p_lsstSwig.i"
%import "lsst/meas/algorithms/algorithmsLib.i"

%lsst_exceptions()

%shared_ptr(lsst::meas::extensions::shapeZ08Katayama::ShapeZ08KatayamaControl);
%shared_ptr(lsst::meas::extensions::shapeZ08Katayama::ShapeZ08KatayamaAlgorithm);
%include "lsst/meas/extensions/shapeZ08Katayama.h"
