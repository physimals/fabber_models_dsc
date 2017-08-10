#pragma once

/*  fwdmodel_dsc.h - DSC model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

/**
 * Base class for the (currently) two DSC models
 *
 * Derived models should implement GetResidual to calculate the residual
 * function using whatever method and parameters they choose. They must
 * also implement GetOptions, Initialize, and GetParameterDefaults to handle their 
 * specific parameters. For these methods, derived models must call the base
 * class methods first to set up common parameters and options (e.g. sig0, cbf
 * parameters and optional dispersion, arterial component, etc..)
 */
class DSCFwdModelBase : public FwdModel
{
public:
    virtual void Initialize(FabberRunData &args);
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key="") const;
    virtual void GetOutputs(std::vector<std::string> &outputs) const;
    
protected:
    /**
     * Calculate the residual function
     *
     * This is implemented separately in the two DSC models
     */
    virtual NEWMAT::ColumnVector CalculateResidual(const NEWMAT::ColumnVector &params) const = 0;

    /**
     * Do a convolution with the AIF.
     *
     * This uses the AIF to form a convolution matrix A and performs m_delt * A * v
     *
     * @param v data to convolute
     * @param aif AIF signal
     */
    NEWMAT::ColumnVector DoConvolution(const NEWMAT::ColumnVector &v, const NEWMAT::ColumnVector &aif) const;
    
    /**
     * Get the AIF for the current voxel - either from suppdata or from user-specified file
     */
    NEWMAT::ColumnVector GetAIF() const;

    /**
     * Shift a time-dependent signal by a fixed delay using linear interpolation
     *
     * Delay is taken in the sense that a positive delay will cause the output
     * signal to be shifted to the right relative to the input signal (i.e. the
     * arrival of the signal is delayed by this amount)
     *
     * For extrapolation, the specified pre-zero value is used if provided, otherwise
     * the input signal is assumed to be initially zero. The signal is assumed to be
     * constant after the final time point.
     *
     * @param sig the time dependent input signal
     * @param delt the time separation between points
     * @param delay the delay to apply (units consistent with delt)
     * @param initial_value Value of the signal before the first time point (default: 0)
     *
     * @return the shifted signal at the same time points as the input signal
     */
    NEWMAT::ColumnVector ApplyDelay(const NEWMAT::ColumnVector &sig, const double dt, const double delay, const double initial_value=0) const;

    /**
     * Apply dispersion of AIF by convolution with a gamma-function based VTF
     */
    NEWMAT::ColumnVector ApplyDispersion(const NEWMAT::ColumnVector &aif, double disp_s, double disp_p) const;

    // scan parameters
    double m_te;
    double m_delt;

    // AIF
    bool m_aifsig;
    NEWMAT::ColumnVector m_aif;

    // Options
    bool m_inferdelay;
    bool m_disp;
    bool m_inferart;
    bool m_art_add_signal;
    std::string m_convmtx;

    virtual void GetParameterDefaults(std::vector<Parameter> &params) const;

    // Lookup the starting indices of the parameters
    // Note that there are two dispersion and two arterial parameters
    int sig0_index() const { return 1; }
    int cbf_index() const { return 2; }
    int delay_index() const { return cbf_index() + (m_inferdelay ? 1 : 0); }
    int disp_index() const { return delay_index() + (m_disp ? 1 : 0); }
    int art_index() const { return disp_index() + (m_disp ? 1 : 0) + (m_inferart ? 1 : 0); }
};

class DSCFwdModel : public DSCFwdModelBase
{
public:
    static FwdModel *NewInstance();

    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;
    
    void Initialize(FabberRunData &args);

protected:
    NEWMAT::ColumnVector CalculateResidual(const NEWMAT::ColumnVector &params) const;

    void GetParameterDefaults(std::vector<Parameter> &params) const;

    int gmu_index() const { return art_index() + (m_inferart ? 1 : 0) + (m_infermtt ? 1 : 0); }
    int lambda_index() const { return gmu_index() + (m_inferlambda ? 1 : 0); }
    int ret_index() const { return lambda_index() + (m_inferret ? 1 : 0); }
    int cbv_index() const { return ret_index() + (m_usecbv ? 1 : 0); }
    
    bool m_infermtt;
    bool m_usecbv;
    bool m_inferlambda;
    bool m_inferret;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DSCFwdModel> registration;
};
