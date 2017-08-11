#!/bin/env python

import os, sys
import traceback

FSLDIR = os.environ["FSLDIR"]
sys.path.insert(0, FSLDIR + "/lib/python")
from fabber import FabberLib, self_test, FabberException

TEST_DATA = {
    "dsc_cpi" : [
        {"te" : 0.065, "aif" : "aif_sig.mat", "delt" : 1.5, "inferdelay" : "", "aifsig" : "", "num-cps" : 5, },
        {"sig0" : [500, 700],
         "cbf" : [0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
         "delay" : [-3, 0, 4],
         "amp1" : [0.1], "amp2" : [0.25], "amp3" : [0.5], "amp4" : [0.5], "amp5" : [0.5]},
        {"nt" : 60} # Other options
    ],
    "dsc" : [
        {"te" : 0.065, "aif" : "aif_sig.mat", "delt" : 1.5, "inferdelay" : "", "aifsig" : "", "infermtt" : "", "inferlambda" : ""},
        {"sig0" : [500, 700],
         "cbf" : [0.01, 0.02, 0.05, 0.1, 0.2, 0.5],
         "delay" : [-3, 0, 4],
         "transitm" : [1.5], "lambda" : [2.3]},
        {"nt" : 60} # Other options
    ],
}

save = "--save" in sys.argv
try:
    for model, test_data in TEST_DATA.items():
        rundata, params, kwargs = test_data
        rundata["model"] = model
        res, log = self_test(model, rundata, params, save_input=save, save_output=save, invert=True, noise=100, **kwargs)
except FabberException, e:
    print e.log
    traceback.print_exc()
except:
    traceback.print_exc()

