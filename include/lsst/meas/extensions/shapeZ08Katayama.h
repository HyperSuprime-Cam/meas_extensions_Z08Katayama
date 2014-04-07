// -*- lsst-c++ -*-

// Feel free to add your own copyright notice here, but note that it must be GPL3-compatible.

// This header file should #include all public header files in the package.
// I've just put all the code here, since I'm only anticipating that we'll need
// this one header, but if you add more, please create a shapeZ08Katayama
// subdirectory and put the additional header files there.  Then #include them
// in this file.

#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace extensions { namespace shapeZ08Katayama {

class ShapeZ08KatayamaAlgorithm;

// This class contains all configurable options that change the details of how the algorithms is run.
// Each option should be defined using the LSST_CONTROL_FIELD macro, which creates a public data member
// as well as some extra hidden things that allow an interface to be created in Python.
class ShapeZ08KatayamaControl : public algorithms::AlgorithmControl {
public:

    // The configuration object must have a default constructor (no arguments).
    ShapeZ08KatayamaControl();

    /// @copydoc lsst::meas::algorithms::AlgorithmControl::clone()
    PTR(ShapeZ08KatayamaControl) clone() const {
        return boost::static_pointer_cast<ShapeZ08KatayamaControl>(_clone());
    }

    /// @copydoc lsst::meas::algorithms::AlgorithmControl::makeAlgorithm()
    PTR(ShapeZ08KatayamaAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        algorithms::AlgorithmMap const & others = algorithms::AlgorithmMap(),
        bool isForced = false
    ) const;

    // Here's an example configuration option.  You can use int, double, std::string, and bool options.
    // You can also create vectors of these (ask me for help if you'd like to do that).
    LSST_CONTROL_FIELD(
        nGrowFootprint, int,
        "Number of pixels to grow the original footprint when determining which pixels to use"
    );

private:
    // boilerplate
    virtual PTR(algorithms::AlgorithmControl) _clone() const;

    // boilerplate
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata,
        algorithms::AlgorithmMap const & others,
        bool isForced
    ) const;
};

class ShapeZ08KatayamaAlgorithm : public algorithms::Algorithm {
public:

    typedef ShapeZ08KatayamaControl Control;

    Control const & getControl() const {
        return static_cast<Control const &>(algorithms::Algorithm::getControl());
    }

    // This constructor is called by the _makeAlgorithm method of the Control class.
    // We could pass the metadata and others map here as well if needed, but I think
    // all we need in this case is the Schema (which we'll use to define the outputs
    // of the algorithm) and the isForced flag, which we'll just use to signal that
    // it doesn't make sense to run this algorithm in forced-photometry mode.
    ShapeZ08KatayamaAlgorithm(afw::table::Schema & schema, Control const & ctrl, bool isForced);

private:

    // This is the main method that you'll have to implement in the source file.
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(ShapeZ08KatayamaAlgorithm);

    // These Key objects are a special kind of pointer for writing the outputs to the catalog.
    // You should have a Key for every value you save, as well as Flag Keys to provide boolean
    // diagnostics.  The main Flag is used to indicate any kind of failure, and should always
    // be set when anythign goes wrong, but you can set additional flags to help describe what
    // went wrong.  Don't worry about creating many Flags; each one only takes up a single bit
    // in each row of the catalog.  The extra flag key here is used to indicate when the object
    // is too close to the edge
    afw::table::Key<double> _e1Key;
    afw::table::Key<double> _e2Key;
    afw::table::Key<afw::table::Flag> _flagKey;
    afw::table::Key<afw::table::Flag> _edgeFlagKey;
};

inline PTR(ShapeZ08KatayamaAlgorithm) ShapeZ08KatayamaControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others,
    bool isForced
) const {
    return boost::static_pointer_cast<ShapeZ08KatayamaAlgorithm>(
        _makeAlgorithm(schema, metadata, others, isForced)
    );
}


}}}} // namespace lsst::meas::extensions::shapeZ08Katayama
