/*  fwdmodel_dsc.h - DSC model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

class DSCFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    virtual void Initialize(ArgsType &args);
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    virtual std::vector<std::string> GetUsage() const;
    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void DumpParameters(const NEWMAT::ColumnVector &vec,
        const string &indents = "") const;

    virtual void NameParams(std::vector<std::string> &names) const;
    virtual int NumParams() const
    {
        return 2 + (infermtt ? 1 : 0) + (inferlambda ? 1 : 0) + (inferdelay ? 1 : 0) + (inferart ? 2 : 0) + (inferret ? 1 : 0) + (usecbv ? 1 : 0) + (dispoption ? 2 : 0);
    }

    virtual ~DSCFwdModel() { return; }
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

protected:
    NEWMAT::ColumnVector aifshift(const NEWMAT::ColumnVector &aif, const double delta, const double hdelt) const;
    void createconvmtx(NEWMAT::LowerTriangularMatrix &A, const NEWMAT::ColumnVector aifnew) const;

    // Constants

    // Lookup the starting indices of the parameters
    int cbf_index() const { return 1; }
    int gmu_index() const { return 1 + (infermtt ? 1 : 0); }
    int lambda_index() const { return 1 + (infermtt ? 1 : 0) + (inferlambda ? 1 : 0); }
    int delta_index() const { return 1 + (infermtt ? 1 : 0) + (inferlambda ? 1 : 0) + (inferdelay ? 1 : 0); }
    int sig0_index() const { return 2 + (infermtt ? 1 : 0) + (inferlambda ? 1 : 0) + (inferdelay ? 1 : 0); }
    int art_index() const { return sig0_index() + (inferart ? 1 : 0); }
    int ret_index() const { return art_index() + (inferart ? 1 : 0) + (inferret ? 1 : 0); } //NB two arterial parameters
    int cbv_index() const { return ret_index() + (usecbv ? 1 : 0); }
    int disp_index() const { return cbv_index() + (dispoption ? 1 : 0); }
    //for ARD
    vector<int> ard_index;

    // scan parameters
    double te;
    double r2;
    double delt;
    NEWMAT::ColumnVector artsig;
    NEWMAT::ColumnVector s;

    bool aifconc;

    bool infermtt;
    bool usecbv;
    bool inferlambda;
    bool inferdelay;
    bool inferart;
    bool artoption;
    bool dispoption;
    bool inferret;
    bool doard;

    bool imageprior;

    std::string convmtx;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, DSCFwdModel> registration;
};
