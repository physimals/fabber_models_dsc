/*  fwdmodel_dsc.cc - Implements a convolution based model for DSC analysis

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dsc.h"

#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"
#include "miscmaths/miscprob.h"

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DSCFwdModel>
    DSCFwdModel::registration("dsc");

static OptionSpec OPTIONS[] = {
    { "te", OPT_FLOAT, "TE echo time in s", OPT_REQ, "" },
    { "delt", OPT_FLOAT, "Time separation between volumes in minutes", OPT_REQ, "0" },
    { "expools", OPT_MATRIX, "ASCII matrix containing extra pool specification", OPT_NONREQ, "" },
    { "ptrain", OPT_MATRIX, "ASCII matrix containing pulsed saturation specification", OPT_NONREQ, "" },
    { "infermtt", OPT_BOOL, "Infer MTT parameter", OPT_NONREQ, "" },
    { "usecbv", OPT_BOOL, "Use CBV", OPT_NONREQ, "" },
    { "inferlambda", OPT_BOOL, "Infer lambda parameter", OPT_NONREQ, "" },
    { "inferdelay", OPT_BOOL, "Infer delay parameter", OPT_NONREQ, "" },
    { "inferart", OPT_BOOL, "Infer arterial component", OPT_NONREQ, "" },
    { "inferret", OPT_BOOL, "Infer RET parameter", OPT_NONREQ, "" },
    { "artoption", OPT_BOOL, "Add signals rather than concentrations", OPT_NONREQ, "" },
    { "disp", OPT_BOOL, "Include some dispersion", OPT_NONREQ, "" },
    { "convmtx", OPT_STR, "Type of convolution matrix: simple or voltera", OPT_NONREQ, "simple" },
    { "aif", OPT_MATRIX, "ASCII matrix containing the arterial signal", OPT_NONREQ, "none" },
    { "aifconc", OPT_BOOL, "Indicates that the AIF is a CTC not signal curve", OPT_NONREQ, "" },
    { "imageprior", OPT_BOOL, "Temp way to indicate we have some image priors (very fixed meaning!)    ", OPT_NONREQ, "" },
    { "" },
};

void DSCFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
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

vector<string> DSCFwdModel::GetUsage() const
{
    vector<string> usage;
    usage.push_back("\nUsage info for --model=dsc:\n");
    usage.push_back("Undefined\n");

    return usage;
}

void DSCFwdModel::Initialize(FabberRunData &args)
{
    // specify command line parameters here
    te = args.GetDouble("te", 0.0);
    delt = args.GetDouble("delt");

    // specify options of the model
    infermtt = args.GetBool("infermtt");
    usecbv = args.GetBool("usecbv");
    if (infermtt & usecbv)
    {
        throw InvalidOptionValue("usecbv, infermtt", "", "Cannot infermtt and usecbv simultaneously");
    }

    inferlambda = args.ReadBool("inferlambda");
    inferdelay = args.ReadBool("inferdelay");
    inferart = args.ReadBool("inferart");   //infer arterial component
    artoption = args.ReadBool("artoption"); //determines if we add concentrations (false) or signals
    inferret = args.ReadBool("inferret");
    dispoption = args.ReadBool("disp"); // determines if we include some dispersion

    convmtx = args.ReadWithDefault("convmtx", "simple");

    // Read in the arterial signal (this will override an image supplied as supplementary data)
    string artfile = args.ReadWithDefault("aif", "none");
    if (artfile != "none")
    {
        artsig = read_ascii_matrix(artfile);
    }

    aifconc = args.ReadBool("aifconc"); // indicates that the AIF is a CTC not signal curve

    doard = false;
    if (inferart)
        doard = true;

    imageprior = args.ReadBool("imageprior"); //temp way to indicate we have some image priors (very fixed meaning!)
}

void DSCFwdModel::DumpParameters(const ColumnVector &vec,
    const string &indent) const
{
}

void DSCFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("cbf");
    if (infermtt)
        names.push_back("transitm");
    if (inferlambda)
        names.push_back("lambda");

    if (inferdelay)
        names.push_back("delay");

    names.push_back("sig0");

    if (inferart)
    {
        names.push_back("abv");
        names.push_back("artdelay");
    }
    if (inferret)
    {
        names.push_back("ret");
    }
    if (usecbv)
    {
        names.push_back("cbv");
    }
    if (dispoption)
    {
        names.push_back("disp_s");
        names.push_back("disp_p");
    }
}

void DSCFwdModel::HardcodedInitialDists(MVNDist &prior,
    MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // CBF
    prior.means(cbf_index()) = 0;
    precisions(cbf_index(), cbf_index()) = 1e-12;
    if (imageprior)
        precisions(cbf_index(), cbf_index()) = 10;

    if (infermtt)
    {
        // Transit mean parameter
        prior.means(gmu_index()) = 1.5;
        precisions(gmu_index(), gmu_index()) = 10;
        if (imageprior)
            precisions(gmu_index(), gmu_index()) = 100;
    }

    if (inferlambda)
    {
        // Transit labmda parameter (log)
        prior.means(lambda_index()) = 2.3;
        precisions(lambda_index(), lambda_index()) = 1;
    }

    if (inferdelay)
    {
        // delay parameter
        prior.means(delta_index()) = 0;
        precisions(delta_index(), delta_index()) = 0.04; //[0.1]; //<1>;
    }

    // signal magnitude parameter
    prior.means(sig0_index()) = 100;
    precisions(sig0_index(), sig0_index()) = 1e-6;
    if (imageprior)
        precisions(sig0_index(), sig0_index()) = 1e12;

    if (inferart)
    {
        //arterial component parameters
        prior.means(art_index()) = 0;
        precisions(art_index(), art_index()) = 1e-12;
        prior.means(art_index() + 1) = 0;
        precisions(art_index() + 1, art_index() + 1) = 0.04;
    }

    if (inferret)
    {
        //some tracer is retained
        prior.means(ret_index()) = 0;
        precisions(ret_index(), ret_index()) = 1e4;
    }

    if (usecbv)
    {
        // CBV is input as an image prior
        prior.means(cbv_index()) = 0;
        precisions(cbv_index(), cbv_index()) = 1e12;
    }

    if (dispoption)
    {
        //including dispersion
        prior.means(disp_index()) = 0.7;
        prior.means(disp_index() + 1) = 0.1;
        precisions(disp_index(), disp_index()) = 100;
        precisions(disp_index() + 1, disp_index() + 1) = 100;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    // Tissue perfusion
    posterior.means(cbf_index()) = 0.1;
    precisions(cbf_index(), cbf_index()) = 0.1;

    if (inferart)
    {
        posterior.means(art_index()) = 0;
        precisions(art_index(), art_index()) = 0.1;
    }

    posterior.SetPrecisions(precisions);
}

void DSCFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
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

    // parameters that are inferred - extract and give sensible names
    float cbf;
    float gmu;    //mean of the transit time distribution
    float lambda; // log lambda (fromt ransit distirbution)
    float cbv;
    float delta;
    float sig0; //'inital' value of the signal
    float artmag;
    float artdelay;
    float tracerret;
    float disp_s;
    float disp_p;

    // extract values from params
    cbf = paramcpy(cbf_index());
    //cbf = exp(params(cbf_index()));
    //if (cbf>1e4) cbf=1e4;

    if (infermtt)
    {
        gmu = params(gmu_index()); //this is the log of the mtt so we can have -ve values
    }
    else
    {
        gmu = 0;
    }

    if (inferlambda)
    {
        lambda = params(lambda_index()); //this is the log of lambda so we can have -ve values
    }
    else
    {
        lambda = 0;
    }

    if (usecbv)
    {
        cbv = params(cbv_index());
    }
    else
    {
        cbv = 0;
    }

    if (inferdelay)
    {
        delta = params(delta_index()); // NOTE: delta is allowed to be negative
    }
    else
    {
        delta = 0;
    }
    sig0 = paramcpy(sig0_index());

    if (inferart)
    {
        artmag = paramcpy(art_index());
        artdelay = params(art_index() + 1);
    }

    if (inferret)
    {
        tracerret = tanh(paramcpy(ret_index()));
    }
    else
        tracerret = 0.0;

    if (dispoption)
    {
        disp_s = params(disp_index());
        disp_p = params(disp_index() + 1);
    }
    else
    {
        disp_s = 0;
        disp_p = 0;
    }

    ColumnVector artsighere; // the arterial signal to use for the analysis
    if (artsig.Nrows() > 0)
    {
        artsighere = artsig; //use the artsig that was loaded by the model
    }
    else
    {
        //use an artsig from supplementary data
        if (suppdata.Nrows() > 0)
        {
            artsighere = suppdata;
        }
        else
        {
            throw FabberRunDataError("No valid AIF provided - require aif option or suppdata");
        }
    }
    // use length of the aif to determine the number of time points
    int ntpts = artsighere.Nrows();

    // sensible limits on delta (beyond which it gets silly trying to estimate it)
    if (delta > ntpts / 2 * delt)
    {
        delta = ntpts / 2 * delt;
    }
    if (delta < -ntpts / 2 * delt)
    {
        delta = -ntpts / 2 * delt;
    }

    //upsampled timeseries
    int upsample;
    int nhtpts;
    float hdelt;
    ColumnVector htsamp;

    // Create vector of sampled times
    ColumnVector tsamp(ntpts);
    for (int i = 1; i <= ntpts; i++)
    {
        tsamp(i) = (i - 1) * delt;
    }

    upsample = 1;
    nhtpts = (ntpts - 1) * upsample + 1;
    htsamp.ReSize(nhtpts);
    htsamp(1) = tsamp(1);
    hdelt = delt / upsample;
    for (int i = 2; i <= nhtpts - 1; i++)
    {
        htsamp(i) = htsamp(i - 1) + hdelt;
    }
    htsamp(nhtpts) = tsamp(ntpts);

    // calculate the arterial input function (from upsampled artsig)
    ColumnVector aif_low(ntpts);
    if (!aifconc)
    {
        aif_low = -1 / te * log(artsighere / artsighere(1)); //using first value from aif input as time zero value
    }
    else
    {
        aif_low = artsighere;
    }

    // upsample the signal
    ColumnVector aif;
    aif.ReSize(nhtpts);
    aif(1) = aif_low(1);
    int j = 0;
    int ii = 0;
    for (int i = 2; i <= nhtpts - 1; i++)
    {
        j = floor((i - 1) / upsample) + 1;
        ii = i - upsample * (j - 1) - 1;
        aif(i) = aif_low(j) + ii / upsample * (aif_low(j + 1) - aif_low(j));
    }
    aif(nhtpts) = aif_low(ntpts);

    // create the AIF matrix - empty for the time being
    LowerTriangularMatrix A(nhtpts);
    A = 0.0;

    // deal with delay parameter - this shifts the aif
    ColumnVector aifnew(aif);
    aifnew = aifshift(aif, delta, hdelt);

    ColumnVector C_art(aif);
    if (inferart)
    {
        //local arterial contribution is the aif, but with a local time shift
        C_art = artmag * aifshift(aif, artdelay, hdelt);
    }

    // Do dispersion of AIF - do this be convolution with a VTF
    if (dispoption)
    {
        ColumnVector vtf;
        vtf.ReSize(nhtpts);
        // Use a gamma VTF
        double s = exp(disp_s);
        double p = exp(disp_p);
        for (int i = 1; i <= nhtpts; i++)
        {
            vtf(i) = pow(s, 1 + s * p) / MISCMATHS::gamma(1 + s * p) * pow(tsamp(i), s * p) * exp(-s * tsamp(i));
        }
        // populate AIF matrix
        createconvmtx(A, aifnew);
        //do the convolution (multiplication)
        aifnew = hdelt * A * vtf;
    }

    // --- Redisue Function ----
    ColumnVector residue;
    residue.ReSize(nhtpts);

    if (lambda > 10)
        lambda = 10;
    if (lambda < -10)
        lambda = -10;
    lambda = exp(lambda);
    if (lambda > 100)
        lambda = 100; //this was 10?

    if (usecbv)
    {
        if (inferart)
        {
            //in this case the cbv image is the 'total' cbv of which part will be arterial
            cbv = cbv - artmag;
        }

        gmu = cbv / (lambda * cbf);
        if (gmu < 1e-6)
            gmu = 1e-6;
    }
    else
    {
        if (gmu > 10)
            gmu = 10;
        if (gmu < -10)
            gmu = -10;
        gmu = exp(gmu);
    }

    float gvar = gmu * gmu / lambda;

    residue = 1 - gammacdf(htsamp.t() - htsamp(1), gmu, gvar).t();
    residue(1) = 1; //always tru - avoid any roundoff errors

    //tracer retention
    residue = (1 - tracerret) * residue + tracerret;

    createconvmtx(A, aifnew);

    // do the multiplication
    ColumnVector C;
    C = cbf * hdelt * A * residue;
    //convert to the DSC signal

    ColumnVector C_low(ntpts);
    for (int i = 1; i <= ntpts; i++)
    {
        C_low(i) = C((i - 1) * upsample + 1);
        if (inferart && !artoption)
        { //add in arterial contribution
            C_low(i) += C_art((i - 1) * upsample + 1);
        }
    }

    ColumnVector sig_art(ntpts);
    result.ReSize(ntpts);
    for (int i = 1; i <= ntpts; i++)
    {
        if (inferart && artoption)
        {
            sig_art(i) = C_art((i - 1) * upsample + 1);
            sig_art(i) = exp(-sig_art(i) * te);
            result(i) = sig0 * (1 + (sig_art(i) - 1) + (exp(-C_low(i) * te) - 1));
        }
        else
        {
            result(i) = sig0 * exp(-C_low(i) * te);
        }
    }

    for (int i = 1; i <= ntpts; i++)
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

FwdModel *DSCFwdModel::NewInstance()
{
    return new DSCFwdModel();
}

ColumnVector DSCFwdModel::aifshift(const ColumnVector &aif, const float delta, const float hdelt) const
{
    // Shift a vector in time by interpolation (linear)
    // NB Makes assumptions where extrapolation is called for.
    int nshift = floor(delta / hdelt);         // number of time points of shift associated with delta
    float minorshift = delta - nshift * hdelt; // shift within the sampled time points (this is always a 'forward' shift)

    ColumnVector aifnew(aif);
    int index;
    int nhtpts = aif.Nrows();
    for (int i = 1; i <= nhtpts; i++)
    {
        index = i - nshift;
        if (index == 1)
        {
            aifnew(i) = aif(1) * minorshift / hdelt;
        } //linear interpolation with zero as 'previous' time point
        else if (index < 1)
        {
            aifnew(i) = 0;
        } // extrapolation before the first time point - assume aif is zero
        else if (index > nhtpts)
        {
            aifnew(i) = aif(nhtpts);
        } // extrapolation beyond the final time point - assume aif takes the value of the final time point
        else
        {
            //linear interpolation
            aifnew(i) = aif(index) + (aif(index - 1) - aif(index)) * minorshift / hdelt;
        }
    }
    return aifnew;
}

void DSCFwdModel::createconvmtx(LowerTriangularMatrix &A, const ColumnVector aifnew) const
{
    // create the convolution matrix
    int nhtpts = aifnew.Nrows();

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
