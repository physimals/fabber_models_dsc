/* dualecho_models.cc Shared library functions for dualecho models

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "fwdmodel_dsc.h"

extern "C" {
	int get_num_models()
	{
		return 1;
	}

	const char *get_model_name(int index)
	{
		switch (index) {
		case 0:
			return "dsc";
			break;
		default:
			return NULL;
		}
	}

	NewInstanceFptr get_new_instance_func(const char *name)
	{
		if (string(name) == "dsc") {
			return DSCFwdModel::NewInstance;
		}
		else {
			return NULL;
		}
	}
}

