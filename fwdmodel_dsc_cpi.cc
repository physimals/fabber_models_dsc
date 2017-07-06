/**
 *  fwdmodel_dsc_cpi.cc
 *
 *  Copyright (C) 2008 University of Oxford
 */

/*  CCOPYRIGHT */

#include "fwdmodel_dsc_cpi.h"
#include "spline_interpolator.h"

#include <fabber_core/easylog.h>

#include "newimage/newimageall.h"
#include <miscmaths/miscprob.h>
#include <newmatio.h>

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace NEWIMAGE;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DSCCpiFwdModel>
    DSCCpiFwdModel::registration("dsc_cpi");

FwdModel *DSCCpiFwdModel::NewInstance()
{
    return new DSCCpiFwdModel();
}

static OptionSpec OPTIONS[] = {
    { "te", OPT_FLOAT, "TE echo time in s", OPT_REQ, "" },
    { "num-cps", OPT_FLOAT, "Number of control points", OPT_REQ, "" },
    { "dt", OPT_FLOAT, "Time separation between volumes in minutes", OPT_REQ, "" },
    { "k", OPT_FLOAT, "CBF scaling parameter", OPT_NONREQ, "1" },
    { "disp", OPT_BOOL, "Apply dispersion to AIF", OPT_NONREQ, "" },
    { "aif", OPT_MATRIX, "ASCII matrix containing the arterial signal", OPT_NONREQ, "none" },
    { "" },
};

void DSCCpiFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
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
    // specify command line parameters here
    m_te = args.GetDouble("te", 0.0);
    m_k = args.GetDouble("k", 0.0);
    m_dt = args.GetDouble("dt");
    m_num_cps = args.GetInt("num_cps", 0);
    m_disp = args.GetBool("disp");
    m_aifsig = args.GetBool("aif-sig");

    // Read in the arterial signal (this will override an image supplied as supplementary data)
    string artfile = args.ReadWithDefault("aif", "none");
    if (artfile != "none")
    {
        m_aif = read_ascii_matrix(artfile);
    }
}

void DSCCpiFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("cbf");
    names.push_back("delay");
    names.push_back("sig0");
    if (m_disp) {
        names.push_back("disp_s");
        names.push_back("disp_p");
    }

    // Amplitude and dt of control points
    for (int i=0; i<m_num_cps; i++) {
        names.push_back("amp" + stringify(i+1));
        names.push_back("dt" + stringify(i+1));
    }
}

int DSCCpiFwdModel::NumParams() const
{
    return 3 + 2*m_num_cps + m_disp ? 2 : 0;
}

void DSCCpiFwdModel::HardcodedInitialDists(MVNDist &prior,
    MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    // Default is mean 0, low precision
    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;
    prior.means = 0;

    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior.SetPrecisions(precisions);
}

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
 * @param dt the time separation between points
 * @param delay the delay to apply (units consistent with dt)
 * @param initial_value Value of the signal before the first time point (default: 0)
 *
 * @return the shifted signal at the same time points as the input signal
 */
ColumnVector DSCCpiFwdModel::apply_delay(const ColumnVector &sig, const double dt, const double delay, const double initial_value) const
{
    // Shift a vector in time by interpolation (linear)
    // NB Makes assumptions where extrapolation is called for.

    // Number of whole time points of shift.
    int nshift = floor(delay / m_dt);  

    // Fractional part of the shift
    double fshift = (delay / m_dt) - nshift; 
    
    ColumnVector signew(sig);
    int index;
    for (int i = 1; i <= data.Nrows(); i++)
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
            // Assume sig is zero before zeroth time point
            signew(i) = initial_value;
        }
        else if (index > data.Nrows())
        {
            // Beyond the final time point - assume signal takes the value of the final time point
            signew(i) = sig(data.Nrows());
        } 
        else
        {
            // Linear interpolation
            signew(i) = sig(index) + (sig(index - 1) - sig(index)) * fshift;
        }
    }
    return signew;
}

void DSCCpiFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // ensure that values are reasonable
    // negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= NumParams(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }

    int nt = data.Nrows();
    double cbf = paramcpy(1); // Must be > 0
    double delay = params(2); // Can be < 0
    double sig0 = params(3); // Can be < 0
    cerr << "CFB: " << cbf << endl;
    cerr << "delay: " << delay << endl;
    cerr << "sig0: " << sig0 << endl;
    int p=4;
    double disp_s, disp_p;
    if (m_disp) {
        double disp_s = params(4);
        double disp_p = params(5);
        p += 2;
        cerr << "disp: s=" << disp_s << ", p=" << disp_p << endl;
    }

    // Control point parameters
    vector<double> cp_ampf, cp_dt;
    for (int i=0; i<m_num_cps; i++) {
        cp_ampf.push_back(params(i*2+p));
        cp_dt.push_back(params(i*2+p+1));
    }

    // Get AIF and apply delay parameter
    ColumnVector aif;
    if (m_aif.Nrows() > 0)
    {
        // AIF was supplied as a parameter
        aif = apply_delay(m_aif, m_dt, delay);
    }
    else if (suppdata.Nrows() > 0)
    {
        // AIF was supplied as suppdata
        aif = apply_delay(suppdata, m_dt, delay);
    }
    else
    {
        throw FabberRunDataError("No valid AIF provided - require aif option or suppdata");
    }
    cerr << "Shifted AIF=" << aif.t();

    if (m_aifsig)
    {
        // AIF is a signal not a concentration 
        aif = -1 / m_te * log(aif / aif(1));
        cerr << "Concentration AIF=" << aif.t();
    }

    if (m_disp)
    {
        // Do dispersion of AIF by convolution with a gamma-function based VTF
        ColumnVector vtf(nt);
        
        double s = exp(disp_s); // FIXME required or use param transforms?
        double p = exp(disp_p);
        for (int i = 1; i <= nt; i++)
        {
            double t = (i - 1) * m_dt;
            cerr << s << ", " << p << endl;
            cerr << pow(s, 1 + s * p) * pow(t, s * p) * exp(-s * t) << endl;
            cerr << MISCMATHS::gamma(2) << endl;
            vtf(i) = pow(s, 1 + s * p) * pow(t, s * p) * exp(-s * t) / MISCMATHS::gamma(1 + s * p);
        }
        cerr << "VTF=" << vtf.t();

        // Do the convolution of AIF with the VTF
        LowerTriangularMatrix A(nt);
        createconvmtx(A, aif);
        aif = m_dt * A * vtf;
        cerr << "Convoluted AIF=" << aif.t();
    }

    // Residual Function is spline-interpolated using the control points
    vector<double> cp_time, cp_amp;
    cp_time.push_back(0);
    cp_amp.push_back(1);
    cerr << "CPs: (0, 1) : ";
    for (int t=0; t<m_num_cps; t++) {
        cp_time.push_back(cp_time[t] + cp_dt[t]);
        cp_amp.push_back(cp_amp[t] * cp_ampf[t]);
        if (t == m_num_cps-1) cp_time[t+1] = (nt-1)*m_dt;
        cerr << "(" << cp_time[t+1] << ", " << cp_amp[t+1] << ") : ";
    }
    cerr << endl;
    
    ColumnVector residual(nt);
    //NaturalSplineInterpolator interp(cp_time, cp_amp);
    PchipInterpolator interp(cp_time, cp_amp);
    for (int i=0; i<nt; i++) {
        residual(i+1) = interp(i*m_dt);
    }
    cerr << "Residue function: " << residual.t();

    // CTC is convolution of AIF with residual
    LowerTriangularMatrix A(nt);
    createconvmtx(A, aif);
    ColumnVector C = cbf * A * residual;
    cerr << "CTC=" << C.t();

    //convert to the DSC signal
    result.ReSize(nt);
    for (int i = 1; i <= nt; i++)
    {
        double sig = sig0 * exp(-C(i) * m_te * m_k);
        result(i) = sig;
    }

    for (int i = 1; i <= nt; i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN of inf in result" << endl;
            LOG << "result: " << result.t() << endl;
            LOG << "params: " << params.t() << endl;
            result = 0.0;
            break;
        }
    }
}

void DSCCpiFwdModel::createconvmtx(LowerTriangularMatrix &A, const ColumnVector aifnew) const
{
    // create the convolution matrix
    int nhtpts = aifnew.Nrows();
    string convmtx = "simple";

    if (convmtx == "simple")
    {
        // Simple convolution matrix
        for (int i = 1; i <= nhtpts; i++)
        {
            for (int j = 1; j <= i; j++)
            {
                A(i, j) = aifnew(i - j + 1); //note we are using the local aifnew here! (i.e. it has been suitably time shifted)
            }
        }
    }
    else if (convmtx == "voltera")
    {
        ColumnVector aifextend(nhtpts + 2);
        ColumnVector zero(1);
        zero = 0;
        aifextend = zero & aifnew & zero;
        int x, y, z;
        //voltera convolution matrix (as defined by Sourbron 2007) - assume zeros outside aif range
        for (int i = 1; i <= nhtpts; i++)
        {
            for (int j = 1; j <= i; j++)
            {
                //cout << i << "  " << j << endl;
                x = i + 1;
                y = j + 1;
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
}
