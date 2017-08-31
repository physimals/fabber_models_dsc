/* dualecho_models.cc Shared library functions for dualecho models

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "dsc_models.h"
#include "fwdmodel_dsc.h"
#include "fwdmodel_dsc_cpi.h"

extern "C" {
int CALL get_num_models()
{
    return 2;
}

const char *CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "dsc";
        break;
    case 1:
        return "dsc_cpi";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "dsc")
    {
        return DSCFwdModel::NewInstance;
    }
    else if (string(name) == "dsc_cpi")
    {
        return DSCCpiFwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
