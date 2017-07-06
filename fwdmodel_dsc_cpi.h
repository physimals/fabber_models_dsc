#pragma once

/** 
 * fwdmodel_dsc_cpi.h - DSC model using control point interpolation
 *
 * Copyright (C) 2007 University of Oxford 
 */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

class DSCCpiFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    virtual void Initialize(ArgsType &args);
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    
    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void NameParams(std::vector<std::string> &names) const;
    virtual int NumParams() const;

    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

protected:

    NEWMAT::ColumnVector apply_delay(const NEWMAT::ColumnVector &sig, const double dt, const double delay, const double initial_value=0) const;
    void createconvmtx(NEWMAT::LowerTriangularMatrix &A, const NEWMAT::ColumnVector aifnew) const;

    // scan parameters
    double m_te;
    double m_k;
    double m_dt;
    int m_num_cps;

    bool m_aifsig;
    bool m_disp;

    NEWMAT::ColumnVector m_aif;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DSCCpiFwdModel> registration;
};

