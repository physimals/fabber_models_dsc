/**
 *  fwdmodel_dsc_cpi.cc
 *
 *  Copyright (C) 2008 University of Oxford
 */

/*  CCOPYRIGHT */

#include "fwdmodel_dsc_cpi.h"
#include "spline_interpolator.h"

#include <fabber_core/easylog.h>
#include <fabber_core/tools.h>
#include <fabber_core/priors.h>
#include <fabber_core/transforms.h>

#include <newimage/newimageall.h>
#include "armawrap/newmat.h"

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DSCCpiFwdModel>
    DSCCpiFwdModel::registration("dsc_cpi");

FwdModel *DSCCpiFwdModel::NewInstance()
{
    return new DSCCpiFwdModel();
}

static OptionSpec OPTIONS[] = {
    { "num-cps", OPT_FLOAT, "Number of control points", OPT_REQ, "" },
    { "infer-cpt", OPT_BOOL, "Infer control point times", OPT_NONREQ, "" },
    { "" },
};

void DSCCpiFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    DSCFwdModelBase::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string DSCCpiFwdModel::GetDescription() const
{
    return "Model for Dynamic Susceptibility Contrast perfusion MRI using Control Point Interpolation method";
}

string DSCCpiFwdModel::ModelVersion() const
{
    string version = "fwdmodel_dsc_cpi.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void DSCCpiFwdModel::Initialize(FabberRunData &args)
{
    DSCFwdModelBase::Initialize(args);

    m_num_cps = args.GetInt("num-cps", 0);
    m_infer_cpt = args.GetBool("infer-cpt");
    m_inc_scaling = args.GetBool("inc-scaling");
}

void DSCCpiFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    DSCFwdModelBase::GetParameterDefaults(params);
    int p=params.size();

    for (int i=0; i<m_num_cps; i++) {
        params.push_back(Parameter(p++, "amp" + stringify(i+1), DistParams(0.5, 1e7), DistParams(0.5, 10), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
        if (m_infer_cpt) {
            params.push_back(Parameter(p++, "cpt" + stringify(i+1), DistParams(1, 10), DistParams(1, 10), PRIOR_NORMAL, TRANSFORM_LOG()));
        }
    }
    if (m_inc_scaling) {
        params.push_back(Parameter(p++, "squash", DistParams(0.5, 1e7), DistParams(0.5, 10), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    }
}

void DSCCpiFwdModel::InitParams(MVNDist &posterior) const
{
    // Voxelwise initialization of posterior.
    assert(posterior.GetSize() == int(m_params.size()));
    int nt = data.Nrows();
    if (nt == 0) return;

    // Set an initial guess for sig0 from the first data point
    posterior.means(sig0_index()) = data(1);

#if 0
    if (m_infer_cpt) {
        // Set initial control points to be evenly spaced
        unsigned int p = cp_index();

        double spacing = double(nt)/m_num_cps;
        for (int i=0; i<m_num_cps-1; i++) {
            //posterior.means(p) = spacing - 0.1;
            //posterior.means(p) = 1.0/m_num_cps - 0.1;
            //posterior.means(p) = 1.0/m_num_cps - 0.1;
            p += 2;
        }
    }
#endif
}

ColumnVector DSCCpiFwdModel::CalculateResidual(const ColumnVector &params) const
{
    int nt = data.Nrows();

    // Control point parameters
    vector<double> cp_time, cp_amp;

    // Always put a control point at the start with amplitude 1
    cp_amp.push_back(1);
    cp_time.push_back(0);

    // The 'amp' parameters are the amplitude reduction factor (<1) at each CP
    // The 'delt' parameters are the time increase from on CP to the next
    int p = cp_index();
    for (int i=0; i<m_num_cps; i++) {
        double ampf = params(p++);
        cp_amp.push_back(cp_amp[i] * ampf);
        //cp_amp.push_back(ampf);
        double tpos = double(i+1)*nt/m_num_cps;
        //#if 0
        if (m_infer_cpt) {
            tpos = cp_time[i] + 0.1 + params(p++);
        }
        //#endif
        #if 0
        // Implementation which allows timepoints to move a bit
        // but not far enough to overlap
        if (m_infer_cpt) {
            tpos += (params(p++) - 0.5) * nt/m_num_cps;
        }
        #endif
        // Global squashing parameter ranging from 0.5 to 2
        if (m_inc_scaling) {
            tpos *= (0.5 + 1.5*params(scaling_index()));
        }
        cp_time.push_back(tpos);
    }

    // Put a control point beyond the end with amplitude zero
    cp_amp.push_back(0);
    double last_time = cp_time[m_num_cps] + 1;
    if (last_time < nt) last_time = nt;
    cp_time.push_back(last_time);

    #if 0
    // Normalize times so they match the data range.
    for (int i=0; i<m_num_cps+1; i++) {
        cp_time[i] = nt*cp_time[i]/cp_time[m_num_cps];
        //cerr << "CP" << i << " - " << cp_amp[i] << ", " << cp_time[i]  << endl;
    }
    #endif

    #if 0
    cerr << "CPs: ";
    for (int i=0; i<cp_amp.size(); i++) {
        cerr << i << " (" << cp_time[i] << ", " << cp_amp[i]  << ") ";
    }
    cerr << endl;
    #endif

    ColumnVector residual(nt);
    //NaturalSplineInterpolator interp(cp_time, cp_amp);
    PchipInterpolator interp(cp_time, cp_amp);
    for (int i=0; i<nt; i++) {
        //cerr << "Interpolating at " << i << endl;
        residual(i+1) = interp(i);
    }
    //cerr << "Residue function: " << residual.t();
    return residual;
}
