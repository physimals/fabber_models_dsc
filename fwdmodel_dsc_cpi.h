#pragma once

/** 
 * fwdmodel_dsc_cpi.h - DSC model using control point interpolation
 *
 * Copyright (C) 2007 University of Oxford 
 */

/*  CCOPYRIGHT */

#include "fwdmodel_dsc.h"

#include <string>
#include <vector>

/**
 * Implementation which calculates the residual using spline interpolation
 *
 * The spline control points are model parameters
 */
class DSCCpiFwdModel : public DSCFwdModelBase
{
public:
    static FwdModel *NewInstance();

    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;
    
    void Initialize(ArgsType &args);
    void InitParams(MVNDist &posterior) const;
    
protected:
    NEWMAT::ColumnVector CalculateResidual(const NEWMAT::ColumnVector &params) const;

    void GetParameterDefaults(std::vector<Parameter> &params) const;

    int cp_index() const { return disp_index() + (m_disp ? 1 : 0) + 1; }
    int scaling_index() const { return cp_index() + (m_infer_cpt ? m_num_cps : 0) + m_num_cps; }

    // Model parameters
    int m_num_cps;
    bool m_infer_cpt;
    bool m_inc_scaling;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DSCCpiFwdModel> registration;
};

