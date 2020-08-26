/*  fwdmodel_dsc.cc - Implements a convolution based model for DSC analysis

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dsc.h"
#include "dscprob/dscprob.h"

#include <fabber_core/easylog.h>
#include <fabber_core/tools.h>
#include <fabber_core/priors.h>

#include <newimage/newimageall.h>
#include <miscmaths/miscprob.h>
#include "armawrap/newmat.h"

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;
using MISCMATHS::gammacdf;

static OptionSpec BASE_OPTIONS[] = {
    { "te", OPT_FLOAT, "TE echo time in s", OPT_REQ, "" },
    { "delt", OPT_FLOAT, "Time separation between volumes in seconds", OPT_REQ, "" },
    { "aif", OPT_MATRIX, "ASCII matrix containing the arterial signal", OPT_NONREQ, "none" },
    { "aifsig", OPT_BOOL, "Indicate that the AIF is a signal curve", OPT_NONREQ, "none" },
    { "aifconc", OPT_BOOL, "Indicate that the AIF is a concentation curve", OPT_NONREQ, "none" },
    { "inferdelay", OPT_BOOL, "Infer delay parameter", OPT_NONREQ, "" },
    { "disp", OPT_BOOL, "Apply dispersion to AIF", OPT_NONREQ, "" },
    { "inferart", OPT_BOOL, "Infer arterial component", OPT_NONREQ, "" },
    { "artoption", OPT_BOOL, "Add signals rather than concentrations", OPT_NONREQ, "" },
    { "convmtx", OPT_STR, "Type of convolution matrix: simple or voltera", OPT_NONREQ, "simple" },
    { "" },
};

void DSCFwdModelBase::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; BASE_OPTIONS[i].name != ""; i++)
    {
        opts.push_back(BASE_OPTIONS[i]);
    }
}

void DSCFwdModelBase::Initialize(FabberRunData &args)
{
    m_te = args.GetDouble("te", 0.0);
    m_delt = args.GetDouble("delt");

    // Read in the arterial signal (this will override an image supplied as supplementary data)
    string artfile = args.ReadWithDefault("aif", "");
    if (artfile != "")
    {
        m_aif = fabber::read_matrix_file(artfile);
    }
    bool aifsig = args.GetBool("aifsig");
    bool aifconc = args.GetBool("aifconc");
    if (aifsig && aifconc) {
        throw InvalidOptionValue("aifconc", "True", "aifsig and aifconc may not be set simultaneously");
    }
    if (!aifsig && !aifconc) {
        WARN_ONCE("WARNING: Neither aifsig nor aifconc set - assuming AIF is a signal curve");
        aifsig = true;
    }
    m_aifsig = aifsig;

    m_disp = args.GetBool("disp");
    m_inferdelay = args.ReadBool("inferdelay");
    m_inferart = args.ReadBool("inferart");
    m_art_add_signal = args.ReadBool("artoption");
    m_convmtx = args.ReadWithDefault("convmtx", "simple");
}

void DSCFwdModelBase::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    params.push_back(Parameter(p++, "sig0", DistParams(100, 1e6), DistParams(100, 1e6)));
    params.push_back(Parameter(p++, "cbf", DistParams(1, 1e6), DistParams(1, 10)));
    if (m_inferdelay) {
        params.push_back(Parameter(p++, "delay", DistParams(0, 25), DistParams(0, 25)));
    }
    if (m_disp) {
        params.push_back(Parameter(p++, "disp_s", DistParams(0, 1e7), DistParams(0, 1e7), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
        params.push_back(Parameter(p++, "disp_p", DistParams(0, 1e7), DistParams(0, 1e7), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    }
    if (m_inferart) {
        params.push_back(Parameter(p++, "abv", DistParams(0, 1e6), DistParams(0, 10)));
        params.push_back(Parameter(p++, "artdelay", DistParams(0, 25), DistParams(0, 25)));
    }
}

void DSCFwdModelBase::GetOutputs(std::vector<std::string> &outputs) const
{
    outputs.push_back("dsc_residual");
}

void DSCFwdModelBase::InitVoxelPosterior(MVNDist &posterior) const
{
    // Initialize signal offset using first data value
    double sig0_data = data(1);
    posterior.means(1) = sig0_data;
}

ColumnVector DSCFwdModelBase::GetAIF() const
{
    ColumnVector aif;
    if (m_aif.Nrows() > 0)
    {
        // AIF was supplied as a parameter
        aif = m_aif;
    }
    else if (suppdata.Nrows() > 0)
    {
        // AIF was supplied as suppdata
        aif = suppdata;
    }
    else
    {
        throw FabberRunDataError("No valid AIF provided - require aif option or suppdata");
    }

    if (aif.Nrows() != data.Nrows()) {
        throw InvalidOptionValue("Length of AIF", stringify(aif.Nrows()), "Expected " + stringify(data.Nrows()));
    }

    if (m_aifsig)
    {
        // AIF is a signal not a concentration
        aif = -1 / m_te * log(aif / aif(1));
    }
    return aif;
}

ColumnVector DSCFwdModelBase::ApplyDelay(const ColumnVector &sig, const double delt, const double delay, const double initial_value) const
{
    // Number of whole time points of shift.
    int nshift = floor(delay / delt);

    // Fractional part of the shift
    double fshift = (delay / delt) - nshift;

    ColumnVector signew(sig);
    int index;
    for (int i = 1; i <= sig.Nrows(); i++)
    {
        index = i - nshift;
        if (index == 1)
        {
            // linear interpolation with initial value as 'previous' time point.
            // Only possible if delay is > 0, so fshift > 0
            signew(i) = sig(1) * (1-fshift) + initial_value*fshift;
        }
        else if (index < 1)
        {
            // Assume sig takes initial_value before zeroth time point
            signew(i) = initial_value;
        }
        else if (index > sig.Nrows())
        {
            // Beyond the final time point - assume signal takes the value of the final time point
            signew(i) = sig(sig.Nrows());
        }
        else
        {
            // Linear interpolation
            signew(i) = sig(index) + (sig(index - 1) - sig(index)) * fshift;
        }
    }
    return signew;
}

// Do dispersion of AIF by convolution with a gamma-function based VTF
ColumnVector DSCFwdModelBase::ApplyDispersion(const ColumnVector &aif, double disp_s, double disp_p) const
{
    unsigned int nt = aif.Nrows();
    ColumnVector vtf(nt);

    double s = exp(disp_s); // FIXME required or use param transforms?
    double p = exp(disp_p);
    for (unsigned int i = 1; i <= nt; i++)
    {
        double t = (i - 1) * m_delt;
        vtf(i) = pow(s, 1 + s * p) * pow(t, s * p) * exp(-s * t) / fabber_dsc::true_gamma(1 + s * p);
    }

    return DoConvolution(vtf, aif);
}

ColumnVector DSCFwdModelBase::DoConvolution(const ColumnVector &v, const ColumnVector &aif) const
{
    // create the convolution matrix
    unsigned int nt = aif.Nrows();
    LowerTriangularMatrix A(nt);
    A = 0;

    if (m_convmtx == "simple")
    {
        // Simple convolution matrix
        for (unsigned int i = 1; i <= nt; i++)
        {
            for (unsigned int j = 1; j <= i; j++)
            {
                A(i, j) = aif(i - j + 1); //note we are using the local aif here! (i.e. it has been suitably time shifted)
            }
        }
    }
    else if (m_convmtx == "voltera")
    {
        ColumnVector aifextend(nt + 2);
        ColumnVector zero(1);
        zero = 0;
        aifextend = zero & aif & zero;
        int x, z;
        //voltera convolution matrix (as defined by Sourbron 2007) - assume zeros outside aif range
        for (unsigned int i = 1; i <= nt; i++)
        {
            for (unsigned int j = 1; j <= i; j++)
            {
                //cout << i << "  " << j << endl;
                x = i + 1;
                z = i - j + 1;
                if (j == 1)
                {
                    A(i, j) = (2 * aifextend(x) + aifextend(x - 1)) / 6;
                }
                else if (j == i)
                {
                    A(i, j) = (2 * aifextend(2) + aifextend(3)) / 6;
                }
                else
                {
                    A(i, j) = (4 * aifextend(z) + aifextend(z - 1) + aifextend(z + 1)) / 6;
                }
            }
        }
    }
    return m_delt * A * v;
}

void DSCFwdModelBase::EvaluateModel(const ColumnVector &params, ColumnVector &result, const std::string &key) const
{
    //cerr << "Params: " << params.t();
    int nt = data.Nrows();
    double sig0 = params(sig0_index());
    double cbf = params(cbf_index());
    double delay = 0;
    if (m_inferdelay) delay = params(delay_index());

    double disp_s=0, disp_p=0;
    if (m_disp) {
        disp_s = params(disp_index());
        disp_p = params(disp_index()+1);
    }

    double artmag=0, artdelay=0;
    if (m_inferart) {
        artmag = params(art_index());
        artdelay = params(art_index()+1);
    }

    // Get original unshifted AIF
    ColumnVector aif_unshifted = GetAIF();

    // Apply delay and dispersion
    ColumnVector aif = ApplyDelay(aif_unshifted, m_delt, delay);
    if (m_disp) aif = ApplyDispersion(aif, disp_s, disp_p);

    ColumnVector C_art(nt);
    if (m_inferart)
    {
        //local arterial contribution is the aif, but with a local time shift
        C_art = artmag * ApplyDelay(aif_unshifted, m_delt, artdelay);
    }

    // Calculation of the residual function
    ColumnVector residual = CalculateResidual(params);

    result.ReSize(nt);
    if (key == "dsc_residual") {
        // Output the residual function
        for (int i = 1; i <= nt; i++)
        {
            result(i) = residual(i);
        }
    }
    else {
        // CTC is convolution of AIF with residual
        ColumnVector C = cbf * DoConvolution(residual, aif);
        //cerr << "CTC=" << C.t();

        //convert to the DSC signal
        for (int i = 1; i <= nt; i++)
        {
            double conc = C(i);
            double sig = 0;
            if (m_inferart) {
                if (m_art_add_signal) {
                    // Add arterial contribution to the signal
                    sig += exp(-C_art(i) * m_te) - 1;
                }
                else {
                    // Add arterial contribution to the concentration
                    conc += C_art(i);
                }
            }

            sig += exp(-conc * m_te) - 1;
            result(i) = sig0 * (1 + sig);
        }
        //cerr << "Result: " << result.t();
    }

    bool warned = false;
    for (int i = 1; i <= nt; i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            if (!warned) {
                LOG << "Warning NaN or inf in result" << endl;
                LOG << "data: " << data.t() << endl;
                LOG << "result: " << result.t() << endl;
                LOG << "params: " << params.t() << endl;
                warned = true;
            }
            result(i) = 0.0;
        }
    }
}

FactoryRegistration<FwdModelFactory, DSCFwdModel>
    DSCFwdModel::registration("dsc");

static OptionSpec OPTIONS[] = {
    { "infermtt", OPT_BOOL, "Infer MTT parameter", OPT_NONREQ, "" },
    { "usecbv", OPT_BOOL, "Use CBV", OPT_NONREQ, "" },
    { "inferlambda", OPT_BOOL, "Infer lambda parameter", OPT_NONREQ, "" },
    { "inferret", OPT_BOOL, "Infer RET parameter", OPT_NONREQ, "" },
    { "" },
};

void DSCFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    DSCFwdModelBase::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string DSCFwdModel::GetDescription() const
{
    return "Model for Dynamic Susceptibility Contrast perfusion MRI ";
}

string DSCFwdModel::ModelVersion() const
{
    string version = "fwdmodel_dsc.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void DSCFwdModel::Initialize(FabberRunData &args)
{
    DSCFwdModelBase::Initialize(args);

    // Options of the model
    m_infermtt = args.GetBool("infermtt");
    m_usecbv = args.GetBool("usecbv");
    if (m_infermtt & m_usecbv)
    {
        throw InvalidOptionValue("usecbv, infermtt", "", "Cannot infermtt and usecbv simultaneously");
    }

    m_inferlambda = args.ReadBool("inferlambda");
    m_inferret = args.ReadBool("inferret");
}

void DSCFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    DSCFwdModelBase::GetParameterDefaults(params);
    int p=params.size();

    if (m_infermtt)
        params.push_back(Parameter(p++, "transitm", DistParams(4.5, 1.1), DistParams(4.5, 1.1), PRIOR_NORMAL, TRANSFORM_LOG()));
    if (m_inferlambda)
        params.push_back(Parameter(p++, "lambda", DistParams(10, 3), DistParams(10, 3), PRIOR_NORMAL, TRANSFORM_LOG()));
    if (m_inferret)
        params.push_back(Parameter(p++, "ret", DistParams(0, 1e-4), DistParams(0, 1e-4)));
    if (m_usecbv)
        params.push_back(Parameter(p++, "cbv", DistParams(0, 1e-6), DistParams(0, 1e-6)));
}

NEWMAT::ColumnVector DSCFwdModel::CalculateResidual(const ColumnVector &params) const
{
    double gmu=0;    // mean of the transit time distribution
    double lambda=0; // lambda (from transit distribution)
    double cbv=0;
    double tracerret=0;

    if (m_infermtt)
    {
        gmu = params(gmu_index());
    }

    if (m_inferlambda)
    {
        lambda = params(lambda_index());
    }

    if (m_usecbv)
    {
        cbv = params(cbv_index());
    }

    if (m_inferret)
    {
        double ret = params(ret_index());
        if (ret < 0) ret = 0;
        tracerret = tanh(ret);
    }

    // Create vector of sampled times for use in dispersion and residue function
    unsigned int nt = data.Nrows();
    ColumnVector tsamp(nt);
    for (unsigned int i = 1; i <= nt; i++)
    {
        tsamp(i) = (i - 1) * m_delt;
    }

    ColumnVector residue(nt);

    if (lambda > 1e5)
        lambda = 1e5;
    if (lambda < 1e-5)
        lambda = 1e-5;

    if (m_usecbv)
    {
        if (m_inferart)
        {
            //in this case the cbv image is the 'total' cbv of which part will be arterial
            cbv = cbv - params(art_index());
        }

        gmu = cbv / (lambda * params(cbf_index()));
        if (gmu < 1e-6)
            gmu = 1e-6;
    }
    else
    {
        if (gmu > 1e5)
            gmu = 1e5;
        if (gmu < 1e-5)
            gmu = 1e-5;
    }

    double gvar = gmu * gmu / lambda;

    residue = 1 - gammacdf(tsamp.t() - tsamp(1), gmu, gvar).t();
    residue(1) = 1; //always true - avoid any roundoff errors

    //tracer retention
    residue = (1 - tracerret) * residue + tracerret;

    return residue;
}

FwdModel *DSCFwdModel::NewInstance()
{
    return new DSCFwdModel();
}
